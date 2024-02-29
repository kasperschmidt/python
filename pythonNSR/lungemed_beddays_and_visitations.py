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
import datetime
from matplotlib.ticker import PercentFormatter

#sys.path.append('C:/Users/kaschm/GitHub/python/pythonNSR/')
import odbc_GetDataFromLPR3 as gdf
import lungemed_beddays_and_visitations as lbv

#-----------------------------------------------------------------------------------------------------------------------
def loadSQL_beddays(filepath='O:\Administration\\02 - Økonomi og PDK\Medarbejdermapper\Kasper\Focus1 - Ad hoc opgaver\Lungemed sengedage og visitationer\Lungemed.sql'):
    """
    Loading SQL query in file "filepath"

    """
    content = open(filepath, 'r').read()
    return content
# -----------------------------------------------------------------------------------------------------------------------
def loadSQL_visitations(filepath='O:\Administration\\02 - Økonomi og PDK\Medarbejdermapper\Kasper\Focus1 - Ad hoc opgaver\Lungemed sengedage og visitationer\Lungemed_visitationsoprindelse_nogroup.sql'):
    """
    Loading SQL query in file "filepath"

    """
    content = open(filepath, 'r').read()
    return content
# -----------------------------------------------------------------------------------------------------------------------
def getdata(verbose=True, filenameext=None):
    """
    Function returning the data structures for the beddays and visitations SQL queries loaded in the load functions.

    Example of yse
    -------
    import lungemed_beddays_and_visitations as lbv
    dataframe_days, dataframe_vis = lbv.getdata()

    import lungemed_beddays_and_visitations as lbv
    dataframe_days, dataframe_vis = lbv.getdata(filenameext='_SUHupdates')
    dataframe_days, dataframe_vis = lbv.getdata(filenameext='_SLAupdates')


    """
    savepath = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Lungemed sengedage og visitationer/plots/'

    if verbose: print(' - Getting LPR3 data from parsing SQL query for "bed days" ')
    savefilename_days  = savepath+'lungemedLPR3_SQLbeddays.xlsx'
    filepath           = savepath+'/../Lungemed.sql'
    if filenameext is not None:
        savefilename_days = savefilename_days.replace('.xlsx', filenameext+'.xlsx')
        filepath          = filepath.replace('.sql', filenameext+'.sql')
    overwrite_days     = True
    dataframe_days     = gdf.returndatapull(lbv.loadSQL_beddays(filepath=filepath), verbose=verbose, savefilename=savefilename_days, overwrite=overwrite_days)

    if verbose: print(' - Getting LPR3 data from parsing SQL query for "visitations" ')
    savefilename_vis   = savepath+'lungemedLPR3_SQLvisitations.xlsx'
    filepath = savepath + '/../Lungemed_visitationsoprindelse_nogroup.sql'
    if filenameext is not None:
        savefilename_vis = savefilename_vis.replace('.xlsx', filenameext+'.xlsx')
        filepath = filepath.replace('.sql', filenameext + '.sql')
    overwrite_vis      = True
    dataframe_vis      = gdf.returndatapull(lbv.loadSQL_visitations(filepath=filepath), verbose=verbose, savefilename=savefilename_vis, overwrite=overwrite_vis)

    return dataframe_days, dataframe_vis

# -----------------------------------------------------------------------------------------------------------------------
def count_occurrences_per_day(measurehours=[8,15,23], untiltoday=False, savedatafile=True, savetype='excel', overwrite=False, verbose=True,
                              filename='lungemedLPR3dataframe', datafileext=None):
    """
    Function to count occurrences per day for various parameters used for bed occupancy and patient statistics

    Example of use
    -------
    import lungemed_beddays_and_visitations as lbv
    df_results = lbv.count_occurrences_per_day(measurehours=[8,15,23], untiltoday=False, savedatafile=True, overwrite=False)

    import lungemed_beddays_and_visitations as lbv
    df_2022results = lbv.count_occurrences_per_day(measurehours=[8,15,23], untiltoday=True, savedatafile=True, filename='lungemedLPR3dataframe_2022data', overwrite=False, datafileext='_updates')

    """
    savepath = 'O:\Administration\\02 - Økonomi og PDK\Medarbejdermapper\Kasper\Focus1 - Ad hoc opgaver\Lungemed sengedage og visitationer\plots\\'
    if os.path.isfile(savepath+filename) and savedatafile and not overwrite:
        sys.exit(' Was asked to store data but overwrite=False and file already exists... hence exiting')

    if verbose: print(' - Getting the data to look at ')
    dataframe_days, dataframe_vis = lbv.getdata(verbose=verbose, filenameext=datafileext)
    outdic = {}

    for measurehour in measurehours:
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if ('SLA' in datafileext) or ('SUH' in datafileext):
            start_day = datetime.datetime.strptime("09-03-2022 " + str(measurehour) + ":00:00", "%d-%m-%Y %H:%M:%S")
        else:
            start_day   = datetime.datetime.strptime("02-02-2019 "+str(measurehour)+":00:00", "%d-%m-%Y %H:%M:%S")

        if untiltoday:
            end_day = datetime.datetime.strptime(str(datetime.datetime.today()).split(' ')[0]+' '+str(measurehour)+":00:00", "%Y-%m-%d %H:%M:%S")
        elif ('SLA' in datafileext) or ('SUH' in datafileext):
            end_day = datetime.datetime.strptime(np.str(dataframe_days['INDTIDSPUNKT_DRGKONTAKT'].max()+datetime.timedelta(days=2)).split(' ')[0]+' '+str(measurehour)+":00:00", "%Y-%m-%d %H:%M:%S")
        else:
            end_day = datetime.datetime.strptime("02-05-2019 " + str(measurehour) + ":00:00", "%d-%m-%Y %H:%M:%S")
        date_list   = [start_day + datetime.timedelta(days=x) for x in range(0, (end_day - start_day).days)]

        if verbose: print(' - Will count how many patients are in beds at any given day between '+
                          start_day.strftime("%d-%m-%Y")+' and '+end_day.strftime("%d-%m-%Y")+' at '+str(measurehour)+" o'clock")

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose: print('---- "Counting parameters from "bed days data frame" ----')
        count_cpr  = [0] * len(date_list)
        occupancy_available  = [0] * len(date_list)
        occupancy_actual     = [0] * len(date_list)

        for pp, patient in enumerate(dataframe_days['CPR']):
            intime  = dataframe_days['INDTIDSPUNKT_DRGKONTAKT'][pp]
            outtime = dataframe_days['UDTIDSPUNKT_DRGKONTAKT'][pp]

            for dd, datecheck in enumerate(np.asarray(date_list)):
                if verbose:
                    infostr = '   Checking the date '+datecheck.strftime("%d-%m-%Y")+' for patient number '+str(pp+1)
                    sys.stdout.write("%s\r" % infostr)
                    sys.stdout.flush()

                if (intime <= datecheck) and (datecheck <= outtime):
                    count_cpr[dd] = count_cpr[dd] + 1

        if verbose: print('\n - Estimating the occupancy in the available and actual beds ')
        for dd, datecheck in enumerate(np.asarray(date_list)):
            if 'SUH' in datafileext:
                if (datecheck > datetime.datetime.strptime("10-03-2022 00:00:00", "%d-%m-%Y %H:%M:%S")) & \
                        (datecheck < datetime.datetime.strptime("27-06-2022 00:00:00", "%d-%m-%Y %H:%M:%S")):
                    NbedsSUH = 18
                elif (datecheck > datetime.datetime.strptime("27-06-2022 00:00:00", "%d-%m-%Y %H:%M:%S")) &\
                        (datecheck < datetime.datetime.strptime("08-08-2022 00:00:00", "%d-%m-%Y %H:%M:%S")): # sommer lavaktivitet
                    NbedsSUH = 14
                elif (datecheck > datetime.datetime.strptime("24-12-2022 00:00:00", "%d-%m-%Y %H:%M:%S")) &\
                        (datecheck < datetime.datetime.strptime("02-01-2023 00:00:00", "%d-%m-%Y %H:%M:%S")): # jul lavaktivitet
                    NbedsSUH = 14
                elif (datecheck > datetime.datetime.strptime("02-01-2023 00:00:00", "%d-%m-%Y %H:%M:%S")) &\
                        (datecheck < datetime.datetime.strptime("01-12-2023 00:00:00", "%d-%m-%Y %H:%M:%S")):
                    NbedsSUH = 20
                elif (datecheck > datetime.datetime.strptime("01-12-2023 00:00:00", "%d-%m-%Y %H:%M:%S")) &\
                        (datecheck < datetime.datetime.strptime("01-12-2024 00:00:00", "%d-%m-%Y %H:%M:%S")):
                    NbedsSUH = 22
                else:
                    NbedsSUH = 18

                occupancy_available[dd] = count_cpr[dd] / NbedsSUH * 100
            elif 'SLA' in datafileext:
                occupancy_available[dd] = count_cpr[dd] / 24. * 100
            else:
                if datecheck < datetime.datetime.strptime("10-06-2021 00:00:00", "%d-%m-%Y %H:%M:%S"):
                    occupancy_available[dd] = count_cpr[dd] / 24. * 100
                else:
                    occupancy_available[dd] = count_cpr[dd] / 16. * 100

                if datecheck < datetime.datetime.strptime("01-03-2021 00:00:00", "%d-%m-%Y %H:%M:%S"):
                    occupancy_actual[dd] = count_cpr[dd] / 24. * 100
                else:
                    occupancy_actual[dd] = count_cpr[dd] / 16. * 100


        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose: print('---- "Counting parameters from "visitations data frame" ----')
        count_vis_aka       = [0] * len(date_list)
        count_vis_lungNAE   = [0] * len(date_list)
        count_vis_lungSLA   = [0] * len(date_list)
        count_vis_other     = [0] * len(date_list)

        for pp, patient in enumerate(dataframe_vis['CPR']):
            intime = dataframe_vis['INDTIDSPUNKT_DRGKONTAKT'][pp]
            outtime = dataframe_vis['UDTIDSPUNKT_DRGKONTAKT'][pp]

            for dd, datecheck in enumerate(np.asarray(date_list)):
                if verbose:
                    infostr = '   Checking the date ' + datecheck.strftime(
                        "%d-%m-%Y") + ' for patient number ' + str(pp + 1)
                    sys.stdout.write("%s\r" % infostr)
                    sys.stdout.flush()

                if intime.strftime("%d-%m-%Y") == datecheck.strftime("%d-%m-%Y"):
                    if ('Akut Afd. 1.sal, Sengeafs., SLA'.lower() in dataframe_vis['SOR_KONTAKT_SP_Afsnit'][pp].lower()) or \
                       ('Akut Afd., Skadestue, SLA'.lower()       in dataframe_vis['SOR_KONTAKT_SP_Afsnit'][pp].lower()) or \
                       ('Akut Afd.stuen, Sengeafs., SLA'.lower()  in dataframe_vis['SOR_KONTAKT_SP_Afsnit'][pp].lower()):
                            count_vis_aka[dd] = count_vis_aka[dd] + 1
                    elif ('Lungemed. Sengeafs., NAE'.lower()  in dataframe_vis['SOR_KONTAKT_SP_Afsnit'][pp].lower()) or \
                         ('Med. Lunge Sengeafs., NAE'.lower() in dataframe_vis['SOR_KONTAKT_SP_Afsnit'][pp].lower()):
                        count_vis_lungNAE[dd] = count_vis_lungNAE[dd] + 1
                    elif ('Med. Lunge Sengeafs., SLA'.lower() in dataframe_vis['SOR_KONTAKT_SP_Afsnit'][pp].lower()):
                        count_vis_lungSLA[dd] = count_vis_lungSLA[dd] + 1
                    else:
                        count_vis_other[dd] = count_vis_other[dd] + 1

        if verbose: print(' - Adding results to output dictionary')
        outdic['dates_'+str(measurehour)]                =  date_list
        outdic['count_cpr_'+str(measurehour)]            =  count_cpr
        outdic['occupancy_available_'+str(measurehour)]  =  occupancy_available
        outdic['occupancy_actual_'+str(measurehour)]     =  occupancy_actual
        outdic['count_vis_aka_'+str(measurehour)]        =  count_vis_aka
        outdic['count_vis_lungNAE_'+str(measurehour)]    =  count_vis_lungNAE
        outdic['count_vis_lungSLA_'+str(measurehour)]    =  count_vis_lungSLA
        outdic['count_vis_other_'+str(measurehour)]      =  count_vis_other

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print('\n - Building data frame and returning count of patients and stats')
    df_results = pd.DataFrame(outdic)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print(' - calculating moving average (30 days window)')
    Ndaysavg = 30
    for measurehour in measurehours:
        df_results['occupancy_available_movingavg_' + str(measurehour)] = \
            df_results['occupancy_available_' + str(measurehour)].rolling(window=Ndaysavg).mean()
        df_results['occupancy_actual_movingavg_' + str(measurehour)] = \
            df_results['occupancy_actual_' + str(measurehour)].rolling(window=Ndaysavg).mean()

    if verbose: print(' - calculating moving average (5 days window)')
    Ndaysavg = 5
    for measurehour in measurehours:
        df_results['occupancy_available_movingavg5days_' + str(measurehour)] = \
            df_results['occupancy_available_' + str(measurehour)].rolling(window=Ndaysavg).mean()
        df_results['occupancy_actual_movingavg5days_' + str(measurehour)] = \
            df_results['occupancy_actual_' + str(measurehour)].rolling(window=Ndaysavg).mean()

    if savedatafile:
        if savetype == 'excel':
            gdf.savefile(df_results, savepath + filename, format='excel', overwrite=overwrite, verbose=verbose)
        else:
            gdf.savefile(df_results, savepath + filename, format='csv', overwrite=overwrite, verbose=verbose)

    return df_results
# -----------------------------------------------------------------------------------------------------------------------
def plot_beddays(loaddatafile='lungemedLPR3_SQLbeddays.xlsx', verbose=True,
                 plotdir='O:\Administration\\02 - Økonomi og PDK\Medarbejdermapper\Kasper\Focus1 - Ad hoc opgaver\Lungemed sengedage og visitationer\plots\\'):
    """
    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    df_results = gdf.loadExcel(plotdir+loaddatafile)
    df_results = df_results.sort_values('INDTIDSPUNKT_DRGKONTAKT')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'beddays_stat.pdf'
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

    xvalues = df_results['INDTIDSPUNKT_DRGKONTAKT']
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d-%Y'))
    plt.gca().xaxis.set_major_locator(mdates.YearLocator())

    xerr = None
    yerr = None

    plt.grid(linestyle=':', linewidth=lthick/2.)
    # --------- RANGES ---------
    ymin = 0.8
    ymax = 150.0
    dy   = ymax - ymin
    plt.ylim([ymin,ymax])

    #plt.xscale('log')
    plt.yscale('log')

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
    pointcolor = cmap(colnorm(40))
    plt.errorbar(xvalues, df_results['KONTAKTDAGE'], xerr=xerr, yerr=yerr,
                 marker='o', lw=0, markersize=marksize, alpha=0.5, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=10.)

    Ndaysavg = 60
    df_results['KONTAKTDAGE_movingavg'] = df_results['KONTAKTDAGE'].rolling(window=Ndaysavg).mean()
    df_results['KONTAKTDAGE_movingavg_std'] = df_results['KONTAKTDAGE'].rolling(window=Ndaysavg).std()

    plt.errorbar(xvalues, df_results['KONTAKTDAGE_movingavg'], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20.,
                 label=str(Ndaysavg)+' dage glidende gennemsnit med standard afvigelse')

    ylow  = df_results['KONTAKTDAGE_movingavg']-df_results['KONTAKTDAGE_movingavg_std']
    yhigh = df_results['KONTAKTDAGE_movingavg']+df_results['KONTAKTDAGE_movingavg_std']
    plt.fill_between(xvalues, ylow, yhigh, color=pointcolor, alpha=0.3, zorder=15)

    lineymin = ymin + dy*0.24
    lineymax = ymin + dy*0.40
    textymin = ymin + dy*0.05
    textymax = ymin + dy*0.10
    plt.plot([datetime.datetime.strptime("10-06-2021", "%d-%m-%Y"), datetime.datetime.strptime("10-06-2021", "%d-%m-%Y")],
             [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
    plt.text(datetime.datetime.strptime("10-06-2021", "%d-%m-%Y"), textymin, '10-06-2021', fontsize=Fsize, rotation=90, color='gray',
             horizontalalignment='center', verticalalignment='bottom')

    plt.plot([datetime.datetime.strptime("01-03-2021", "%d-%m-%Y"), datetime.datetime.strptime("01-03-2021", "%d-%m-%Y")],
             [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
    plt.text(datetime.datetime.strptime("01-03-2021", "%d-%m-%Y"), textymin, '01-03-2021', fontsize=Fsize, rotation=90, color='gray',
             horizontalalignment='center', verticalalignment='bottom')

    # --------- LABELS ---------
    plt.xlabel('Dato for måling')
    plt.ylabel('Kontaktdage Lungemd. NAE')

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=2, numpoints=1,
                     bbox_to_anchor=(0.5, 1.12), )  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------

    plt.savefig(plotdir+plotname)
    plt.clf()
    plt.close('all')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# -----------------------------------------------------------------------------------------------------------------------
def load_occupancy_updates(verbose=True):
    """

    """
    path_updates = 'O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Lungemed sengedage og visitationer/månedligetræk fra SP/'
    updatenames = 'belægning_lungemed_NAE_*.xlsx'
    filelist_updates  = glob.glob(path_updates + updatenames)
    updates_framelist = []
    if verbose: print(' - Loading dataframes with updates from ' + updatenames)
    for file_update in filelist_updates:
        if verbose: print('    - including '+file_update.split('/')[-1])
        df_update = pd.read_excel(file_update, sheet_name='Belægningsprocenter')
        updates_framelist = updates_framelist + [df_update]
    df_updates = pd.concat(updates_framelist, ignore_index=True)
    if verbose: print(' - returning concatenation of '+ str(len(updates_framelist)) +' dataframes for updates')
    return df_updates

# -----------------------------------------------------------------------------------------------------------------------
def plot_perday_occupancy(measurehours=[23], loaddatafile='lungemedLPR3dataframe.xlsx', concatenateLPR=False, verbose=True, untiltoday=True, include_updates=True,
                          plotdir='O:/Administration/02 - Økonomi og PDK/Medarbejdermapper/Kasper/Focus1 - Ad hoc opgaver/Lungemed sengedage og visitationer/plots/',
                          datemin_mostrecent = '01-01-2022'):
    """

    Parameters
    ----------
    measurehours   The hours at which the occupancy should be evaluated
    loaddatafile   If loaddatafile is not None, the script will attempt to load the data from a csv file instead of
                   pulling it from LPR3 via SQL and performing the calculations. The file is assumed to be loaceted in
                   plotdir
    verbose        Toggle verbosity
    untiltoday     Grab data from 02-02-2019 (start of LPR3) until today. If False only a short period after LPR3 start
                   is considered for testing.
    plotdir        Directory to put plots in

    Returns
    -------
    Plots of the pulled data

    Example of use
    -------
    import lungemed_beddays_and_visitations as lbv
    lbv.plot_perday_occupancy(measurehours=[8, 15, 23], untiltoday=True, loaddatafile='lungemedLPR3dataframe.xlsx')

    """
    if verbose: print(' - Generating data for plots')
    if loaddatafile is None:
        df_results = lbv.count_occurrences_per_day(measurehours=measurehours, untiltoday=untiltoday, verbose=verbose, savedatafile=False)
    else:
        if concatenateLPR:
            dffiles = glob.glob(plotdir+loaddatafile.replace('.xlsx','*.xlsx'))
            df_list = []
            for dffile in dffiles:
                if verbose: print(' - Load Excel file with dataframe:\n   '+dffile)
                df = pd.read_excel(dffile)
                df_list.append(df)
            if verbose: print(' - Concatenating the '+str(len(dffiles))+' dataframes found from globbing for \n   '+plotdir+loaddatafile.replace('.xlsx','*.xlsx'))
            df_results = pd.concat(df_list, ignore_index=True)
            #df_results = df_results.sort_values(by='dates_23')
        else:
            df_results = pd.read_excel(plotdir+loaddatafile)

            df_SUHupdates = pd.read_excel(plotdir+loaddatafile.replace('.xlsx','_SUHupdates.xlsx'))
            df_SLAupdates = pd.read_excel(plotdir+loaddatafile.replace('.xlsx','_SLAupdates.xlsx'))

    if verbose: print(' - Removing incomplete data in baseline data frame')
    datemin    = datetime.datetime.strptime('01-11-2021', "%d-%m-%Y")
    datetrans  = datetime.datetime.strptime('09-03-2022', "%d-%m-%Y") # Date for transfer to SUH
    dropval    = np.where((df_results['dates_23'] > datemin) & (df_results['dates_23'] < datetrans))[0]
    df_results = df_results.drop(df_results.index[dropval])

    if include_updates:
        df_updates     = load_occupancy_updates()
        udates         = np.unique(df_updates['Belægning kalenderdato'].values)
        days_updates   = np.sort(udates[~np.isnan(udates)])
        datemin        = datetime.datetime.strptime(datemin_mostrecent, "%d-%m-%Y")
        dropval        = np.where(df_updates['Belægning kalenderdato'] < datemin)[0]
        df_lastmonth   = df_updates.drop(df_updates.index[dropval])
        udates         = np.unique(df_lastmonth['Belægning kalenderdato'].values)
        days_lastmonth = np.sort(udates[~np.isnan(udates)])

        occupancy_updates = {}
        occupancy_updates_movingavg = {}
        occupancy_lastmonth = {}
        for measurehour in measurehours:
            Ndaysavg = 30
            hourstr = str("%.2d" % measurehour)+':00'

            datestrarr = np.asarray([dstr[11:16] for dstr in np.asarray(df_updates['Belægning tidspunkt'].values).astype(str)])
            occupdates = df_updates['Belægning i %'].values[datestrarr == hourstr]
            occupancy_updates['occupancy_updates_' + str(measurehour)] = occupdates[~np.isnan(occupdates)]

            occupancy_updates_movingavg['updates_movingavg_'+str(measurehour)] = pd.DataFrame(occupancy_updates['occupancy_updates_'+str(measurehour)]).rolling(window=Ndaysavg).mean()

            datestrarr = np.asarray([dstr[11:16] for dstr in np.asarray(df_lastmonth['Belægning tidspunkt'].values).astype(str)])
            occupancy_lastmonth['occupancy_lastmonth_' + str(measurehour)]     = df_lastmonth['Belægning i %'].values[datestrarr == hourstr]
            occupancy_lastmonth['occupancy_lastmonth_' + str(measurehour)] = \
            occupancy_lastmonth['occupancy_lastmonth_' + str(measurehour)][
                ~np.isnan(occupancy_lastmonth['occupancy_lastmonth_' + str(measurehour)])]

        occupancy_updates_movingavg5days = {}
        for measurehour in measurehours:
            Ndaysavg = 5
            occupancy_updates_movingavg['updates_movingavg5days_'+str(measurehour)] = pd.DataFrame(occupancy_updates['occupancy_updates_'+str(measurehour)]).rolling(window=Ndaysavg).mean()


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for measurehour in measurehours:
        plotname = 'occupancy_kl'+str(measurehour)+'.pdf'
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

        xvalues = df_results['dates_'+str(measurehour)]
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d-%Y'))
        plt.gca().xaxis.set_major_locator(mdates.YearLocator())
        plt.grid(color='gray', linestyle=':', linewidth=1, alpha=0.5)

        xerr = None
        yerr = None
        # --------- RANGES ---------
        ymin = 10.0
        ymax = 150.0
        dy = ymax - ymin
        plt.ylim([ymin, ymax])

        # --------- COLORMAP ---------
        cmap = plt.cm.get_cmap('viridis')
        #cmap = plt.cm.get_cmap('plasma')
        cmin = 0
        cmax = 100

        colnorm = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax)
        cmaparr = np.linspace(cmin, cmax, num=50)
        m = plt.cm.ScalarMappable(cmap=cmap)
        m.set_array(cmaparr)
        #cb = plt.colorbar(m)
        #cb.set_label('belægningsprocent')

        # --------- POINT AND CURVES ---------

        ##### Baseline data #####
        pointcolor = cmap(colnorm(35))
        lprbaseline = np.where((df_results['dates_23'] < datetrans))[0]
        plt.errorbar(xvalues[df_results.index[lprbaseline]], df_results['occupancy_available_'+str(measurehour)][df_results.index[lprbaseline]], xerr=xerr, yerr=yerr,
                     marker='o', lw=0, markersize=marksize, alpha=0.5,
                     markerfacecolor=pointcolor, ecolor=pointcolor,
                     markeredgecolor=pointcolor, zorder=10,
                     label='NAE baseline (2019-2021)')
        plt.errorbar(xvalues[df_results.index[lprbaseline]], df_results['occupancy_available_movingavg_'+str(measurehour)][df_results.index[lprbaseline]], xerr=xerr, yerr=yerr,
                     marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                     markerfacecolor=pointcolor, ecolor=pointcolor,
                     markeredgecolor=pointcolor, zorder=20.,
                     label='30 dage glidende gennemsnit')

        pointcolor = cmap(colnorm(45))
        datemin_corr = datetime.datetime.strptime('01-03-2021', "%d-%m-%Y")
        datemax_corr = datetime.datetime.strptime('10-06-2021', "%d-%m-%Y")
        corrent      = np.where((df_results['dates_23'] > datemin_corr) & (df_results['dates_23'] < datemax_corr))[0]

        plt.errorbar(xvalues[df_results.index[corrent]], df_results['occupancy_actual_'+str(measurehour)][df_results.index[corrent]], xerr=xerr, yerr=yerr,
                     marker='o', lw=0, markersize=marksize, alpha=0.5,
                     markerfacecolor=pointcolor, ecolor=pointcolor,
                     markeredgecolor=pointcolor, zorder=8.,
                     label='NAE "Reelle" sengekapacitet')
        plt.errorbar(xvalues[df_results.index[lprbaseline]], df_results['occupancy_actual_movingavg_'+str(measurehour)][df_results.index[lprbaseline]], xerr=xerr, yerr=yerr,
                     marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                     markerfacecolor=pointcolor, ecolor=pointcolor,
                     markeredgecolor=pointcolor, zorder=18.) #label='30 dage glidende gennemsnit'

        ##### LPR3 updates #####
        if concatenateLPR:
            pointcolor = cmap(colnorm(55))
            lprupdates      = np.where((df_results['dates_23'] > datetrans))[0]
            plt.errorbar(xvalues[df_results.index[lprupdates]], df_results['occupancy_actual_'+str(measurehour)][df_results.index[lprupdates]], xerr=xerr, yerr=yerr,
                         marker='o', lw=0, markersize=marksize, alpha=0.5,
                         markerfacecolor=pointcolor, ecolor=pointcolor,
                         markeredgecolor=pointcolor, zorder=8.,
                         label='"LPR3 updates" sengekapacitet')

        else:
            pointcolor = cmap(colnorm(65))
            plt.errorbar(df_SLAupdates['dates_'+str(measurehour)], df_SLAupdates['occupancy_available_'+str(measurehour)], xerr=xerr, yerr=yerr,
                         marker='o', lw=0, markersize=marksize, alpha=0.5,
                         markerfacecolor=pointcolor, ecolor=pointcolor,
                         markeredgecolor=pointcolor, zorder=8.,
                         label='SLA LPR3 (siden 10-03-2022)')
            plt.errorbar(df_SLAupdates['dates_'+str(measurehour)], df_SLAupdates['occupancy_available_movingavg_' + str(measurehour)], xerr=xerr, yerr=yerr,
                         marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                         markerfacecolor=pointcolor, ecolor=pointcolor,
                         markeredgecolor=pointcolor, zorder=20.) #label='30 dage glidende gennemsnit'

            pointcolor = cmap(colnorm(85))
            plt.errorbar(df_SUHupdates['dates_'+str(measurehour)], df_SUHupdates['occupancy_available_'+str(measurehour)], xerr=xerr, yerr=yerr,
                         marker='o', lw=0, markersize=marksize, alpha=0.5,
                         markerfacecolor=pointcolor, ecolor=pointcolor,
                         markeredgecolor=pointcolor, zorder=8.,
                         label='SUH LPR3 (siden 10-03-2022)')
            plt.errorbar(df_SUHupdates['dates_'+str(measurehour)], df_SUHupdates['occupancy_available_movingavg_' + str(measurehour)], xerr=xerr, yerr=yerr,
                         marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                         markerfacecolor=pointcolor, ecolor=pointcolor,
                         markeredgecolor=pointcolor, zorder=20.) #label='30 dage glidende gennemsnit'

        ##### SP data #####
        if include_updates:
            pointcolor = 'black'
            plt.errorbar(days_updates, occupancy_updates['occupancy_updates_' + str(measurehour)],
                         xerr=xerr, yerr=yerr, marker='o', lw=0, markersize=marksize, alpha=0.5,
                         markerfacecolor=pointcolor, ecolor=pointcolor,
                         markeredgecolor=pointcolor, zorder=8.,
                         label='NAE 01-12-2021 til 10-03-2022')
            plt.errorbar(days_updates, occupancy_updates_movingavg['updates_movingavg_'+str(measurehour)], xerr=xerr, yerr=yerr,
                         marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                         markerfacecolor=pointcolor, ecolor=pointcolor,
                         markeredgecolor=pointcolor, zorder=18.)#,
                         #label='30 dage glidende gennemsnit')

            pointcolor = 'black'
            plt.errorbar(days_lastmonth, occupancy_lastmonth['occupancy_lastmonth_' + str(measurehour)], xerr=xerr, yerr=yerr,
                         marker='o', lw=0, markersize=marksize, alpha=0.5,
                         markerfacecolor=pointcolor, ecolor=pointcolor,
                         markeredgecolor=pointcolor, zorder=30.) # label='Seneste måned (siden '+datemin_mostrecent+')'

        plt.plot(xvalues, np.zeros(len(xvalues)) + 100, '--', color='black', lw=lthick, zorder=5)
        plt.plot(xvalues, np.zeros(len(xvalues)) + 85, ':', color='black', lw=lthick, zorder=5)

        lineymin = ymin + dy * 0.24
        lineymax = ymin + dy * 0.40
        textymin = ymin + dy * 0.05
        textymax = ymin + dy * 0.10
        plt.plot([datetime.datetime.strptime("10-06-2021", "%d-%m-%Y"),
                  datetime.datetime.strptime("10-06-2021", "%d-%m-%Y")],
                 [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
        plt.text(datetime.datetime.strptime("10-06-2021", "%d-%m-%Y"), textymin, '10-06-2021', fontsize=Fsize,
                 rotation=90, color='gray',
                 horizontalalignment='center', verticalalignment='bottom')

        plt.plot([datetime.datetime.strptime("01-03-2021", "%d-%m-%Y"),
                  datetime.datetime.strptime("01-03-2021", "%d-%m-%Y")],
                 [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
        plt.text(datetime.datetime.strptime("01-03-2021", "%d-%m-%Y"), textymin, '01-03-2021', fontsize=Fsize,
                 rotation=90, color='gray',
                 horizontalalignment='center', verticalalignment='bottom')

        # sommer lavaktivitet
        plt.plot([datetime.datetime.strptime("27-06-2022", "%d-%m-%Y"),
                  datetime.datetime.strptime("27-06-2022", "%d-%m-%Y")],
                 [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
        plt.plot([datetime.datetime.strptime("08-08-2022", "%d-%m-%Y"),
                  datetime.datetime.strptime("08-08-2022", "%d-%m-%Y")],
                 [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
        plt.text(datetime.datetime.strptime("15-07-2022", "%d-%m-%Y"), textymin, 'SUH 14 senge', fontsize=Fsize,
                 rotation=90, color='gray',
                 horizontalalignment='center', verticalalignment='bottom')

        # jul lavaktivitet
        plt.plot([datetime.datetime.strptime("24-12-2022", "%d-%m-%Y"),
                  datetime.datetime.strptime("24-12-2022", "%d-%m-%Y")],
                 [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
        plt.plot([datetime.datetime.strptime("02-01-2023", "%d-%m-%Y"),
                  datetime.datetime.strptime("02-01-2023", "%d-%m-%Y")],
                 [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
        plt.text(datetime.datetime.strptime("29-12-2022", "%d-%m-%Y"), textymin, 'SUH 14 senge', fontsize=Fsize,
                 rotation=90, color='gray',
                 horizontalalignment='center', verticalalignment='bottom')

        # --------- LABELS ---------
        plt.xlabel('Dato for måling (kl '+str(measurehour)+')')
        plt.ylabel('Belægningsprocent Lungemedicin')

        # --------- RANGES ---------
        #xmin = np.min(xvalues[np.isfinite(xvalues)])
        #xmax = np.max(xvalues[np.isfinite(xvalues)])
        #dx = xmax - xmin

        #ymin = np.min(df_results['occupancy_available'][np.isfinite(df_results['occupancy_available'])])
        #ymax = np.max(df_results['occupancy_available'][np.isfinite(df_results['occupancy_available'])])
        #dy = ymax - ymin

        #plt.xlim([xmin - dx * 0.05, xmax + dx * 0.05])
        #plt.ylim([ymin - dy * 0.05, ymax + dy * 0.05])
        #plt.ylim([0,150])

        # if logx:
        #     plt.xscale('log')
        # if logy:
        #     plt.yscale('log')

        plt.ylim(0, plt.gca().get_ylim()[1])
        plt.gca().yaxis.set_major_formatter(PercentFormatter(1*100.))
        # --------- LEGEND ---------
        #plt.errorbar(-5000, -5000, xerr=None, yerr=1, marker='o', lw=0, markersize=marksize, alpha=1.0,
        #             markerfacecolor='k', ecolor='k', markeredgecolor='black', zorder=1, label='MUSE-Wide LAE')
        #plt.errorbar(-5000, -5000, xerr=None, yerr=None, marker='*', lw=0, markersize=marksize * 2, alpha=1.0,
        #             markerfacecolor='None', ecolor='None', markeredgecolor='black', zorder=1, label='AGN')
        #plt.errorbar(-5000, -5000, xerr=None, yerr=None, marker='D', lw=0, markersize=marksize, alpha=1.0,
        #             markerfacecolor='None', ecolor='None', markeredgecolor='black', zorder=1, label='AGN candidate')

        leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=3, numpoints=1,
                         bbox_to_anchor=(0.47, 1.12), )  # add the legend
        leg.get_frame().set_alpha(0.7)
        # --------------------------

        plt.savefig(plotdir+plotname)
        plt.clf()
        plt.close('all')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for measurehour in measurehours:
        plotname = 'occupancy_zoom_kl'+str(measurehour)+'.pdf'
        if verbose: print(' - Initiating '+plotname)

        fig = plt.figure(figsize=(9, 6))
        fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.1, right=0.97, bottom=0.15, top=0.90)
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

        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d-%m-%y'))
        plt.gca().xaxis.set_major_locator(mdates.MonthLocator())
        plt.grid(color='gray', linestyle=':', linewidth=1, alpha=0.5)

        xerr = None
        yerr = None
        # --------- RANGES ---------
        ymin = 10.0
        ymax = 150.0
        dy = ymax - ymin
        plt.ylim([ymin, ymax])

        # --------- COLORMAP ---------
        cmap = plt.cm.get_cmap('viridis')
        #cmap = plt.cm.get_cmap('plasma')
        cmin = 0
        cmax = 100

        colnorm = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax)
        cmaparr = np.linspace(cmin, cmax, num=50)
        m = plt.cm.ScalarMappable(cmap=cmap)
        m.set_array(cmaparr)
        #cb = plt.colorbar(m)
        #cb.set_label('belægningsprocent')

        # --------- POINT AND CURVES ---------
        ##### SP data #####
        showupdates = False
        if showupdates:
            pointcolor = 'black'
            # plt.errorbar(days_updates, occupancy_updates['occupancy_updates_' + str(measurehour)],
            #              xerr=xerr, yerr=yerr, marker='o', lw=0, markersize=marksize, alpha=0.5,
            #              markerfacecolor=pointcolor, ecolor=pointcolor,
            #              markeredgecolor=pointcolor, zorder=8.,
            #              label='NAE 01-12-2021 til 10-03-2022')
            plt.errorbar(days_updates, occupancy_updates_movingavg['updates_movingavg5days_'+str(measurehour)], xerr=xerr, yerr=yerr,
                         marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                         markerfacecolor=pointcolor, ecolor=pointcolor,
                         markeredgecolor=pointcolor, zorder=18.)#,
                         #label='30 dage glidende gennemsnit')

            pointcolor = 'black'
            plt.errorbar(days_lastmonth, occupancy_lastmonth['occupancy_lastmonth_' + str(measurehour)], xerr=xerr, yerr=yerr,
                         marker='o', lw=0, markersize=marksize, alpha=0.5,
                         markerfacecolor=pointcolor, ecolor=pointcolor,
                         markeredgecolor=pointcolor, zorder=30.) # label='Seneste måned (siden '+datemin_mostrecent+')'

        ##### LPR3 updates #####
        pointcolor = cmap(colnorm(65))
        # plt.errorbar(df_SLAupdates['dates_'+str(measurehour)], df_SLAupdates['occupancy_available_'+str(measurehour)], xerr=xerr, yerr=yerr,
        #              marker='o', lw=0, markersize=marksize, alpha=0.5,
        #              markerfacecolor=pointcolor, ecolor=pointcolor,
        #              markeredgecolor=pointcolor, zorder=8.,
        #              label='SLA LPR3 (siden 10-03-2022)')
        plt.errorbar(df_SLAupdates['dates_'+str(measurehour)], df_SLAupdates['occupancy_available_movingavg5days_' + str(measurehour)], xerr=xerr, yerr=yerr,
                     marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                     markerfacecolor=pointcolor, ecolor=pointcolor,
                     markeredgecolor=pointcolor, zorder=20.,
                     label='SLA LPR3 (siden 10-03-2022)') #label='30 dage glidende gennemsnit'

        pointcolor = cmap(colnorm(85))
        # plt.errorbar(df_SUHupdates['dates_'+str(measurehour)], df_SUHupdates['occupancy_available_'+str(measurehour)], xerr=xerr, yerr=yerr,
        #              marker='o', lw=0, markersize=marksize, alpha=0.5,
        #              markerfacecolor=pointcolor, ecolor=pointcolor,
        #              markeredgecolor=pointcolor, zorder=8.,
        #              label='SUH LPR3 (siden 10-03-2022)')
        plt.errorbar(df_SUHupdates['dates_'+str(measurehour)], df_SUHupdates['occupancy_available_movingavg5days_' + str(measurehour)], xerr=xerr, yerr=yerr,
                     marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                     markerfacecolor=pointcolor, ecolor=pointcolor,
                     markeredgecolor=pointcolor, zorder=20.,
                     label='SUH LPR3 (siden 10-03-2022)') #label='30 dage glidende gennemsnit'

        ##### SUH 14 to 18 bed adjustment #####
        pointcolor = 'red'
        ent_adjust18 = np.where((df_SUHupdates['dates_' + str(measurehour)] > datetime.datetime.strptime("27-06-2022 00:00:00", "%d-%m-%Y %H:%M:%S")) & \
                               (df_SUHupdates['dates_' + str(measurehour)] < datetime.datetime.strptime("08-08-2022 00:00:00", "%d-%m-%Y %H:%M:%S")))[0]
        plt.errorbar(df_SUHupdates['dates_' + str(measurehour)][ent_adjust18] ,
                     df_SUHupdates['occupancy_available_movingavg5days_' + str(measurehour)][ent_adjust18]*14./18., xerr=xerr, yerr=yerr,
                     marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                     markerfacecolor=pointcolor, ecolor=pointcolor,
                     markeredgecolor=pointcolor, zorder=20.,
                     label='SUH justeret fra 14 til 18 senge ')

        pointcolor = 'red'
        ent_adjust18 = np.where((df_SUHupdates['dates_' + str(measurehour)] > datetime.datetime.strptime("24-12-2022 00:00:00", "%d-%m-%Y %H:%M:%S")) & \
                               (df_SUHupdates['dates_' + str(measurehour)] < datetime.datetime.strptime("02-01-2023 00:00:00", "%d-%m-%Y %H:%M:%S")))[0]
        plt.errorbar(df_SUHupdates['dates_' + str(measurehour)][ent_adjust18] ,
                     df_SUHupdates['occupancy_available_movingavg5days_' + str(measurehour)][ent_adjust18]*14./18., xerr=xerr, yerr=yerr,
                     marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                     markerfacecolor=pointcolor, ecolor=pointcolor,
                     markeredgecolor=pointcolor, zorder=20.)
        #########################################

        xrangevals = [df_SUHupdates['dates_' + str(measurehour)].values[0], df_SUHupdates['dates_' + str(measurehour)].values[-1]]
        #xrangevals = [days_updates[0], df_SUHupdates['dates_' + str(measurehour)].values[-1]]
        plt.plot(xrangevals,np.zeros(len(xrangevals)) + 100, '--', color='black', lw=lthick, zorder=5)
        plt.plot(xrangevals, np.zeros(len(xrangevals)) + 85, ':', color='black', lw=lthick, zorder=5)

        lineymin = ymin + dy * 0.24
        lineymax = ymin + dy * 0.40
        textymin = ymin + dy * 0.05
        textymax = ymin + dy * 0.10

        # sommer lavaktivitet
        plt.plot([datetime.datetime.strptime("27-06-2022", "%d-%m-%Y"),
                  datetime.datetime.strptime("27-06-2022", "%d-%m-%Y")],
                 [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
        plt.plot([datetime.datetime.strptime("08-08-2022", "%d-%m-%Y"),
                  datetime.datetime.strptime("08-08-2022", "%d-%m-%Y")],
                 [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)

        plt.text(datetime.datetime.strptime("10-07-2022", "%d-%m-%Y"), textymin, 'SUH 14 senge', fontsize=Fsize,
                 rotation=90, color='gray',
                 horizontalalignment='center', verticalalignment='bottom')

        plt.text(datetime.datetime.strptime("23-07-2022", "%d-%m-%Y"), textymin, 'SUH 18 senge', fontsize=Fsize,
                 rotation=90, color='red',
                 horizontalalignment='center', verticalalignment='bottom')


        # jul lavaktivitet
        plt.plot([datetime.datetime.strptime("24-12-2022", "%d-%m-%Y"),
                  datetime.datetime.strptime("24-12-2022", "%d-%m-%Y")],
                 [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
        plt.plot([datetime.datetime.strptime("02-01-2023", "%d-%m-%Y"),
                  datetime.datetime.strptime("02-01-2023", "%d-%m-%Y")],
                 [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)

        plt.text(datetime.datetime.strptime("27-12-2022", "%d-%m-%Y"), textymin, 'SUH 14 senge', fontsize=Fsize,
                 rotation=90, color='gray',
                 horizontalalignment='center', verticalalignment='bottom')

        plt.text(datetime.datetime.strptime("30-12-2022", "%d-%m-%Y"), textymin, 'SUH 18 senge', fontsize=Fsize,
                 rotation=90, color='red',
                 horizontalalignment='center', verticalalignment='bottom')


        # SUH opnormerer til 22 senge
        plt.plot([datetime.datetime.strptime("01-12-2023", "%d-%m-%Y"),
                  datetime.datetime.strptime("01-12-2023", "%d-%m-%Y")],
                 [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)

        plt.text(datetime.datetime.strptime("15-12-2023", "%d-%m-%Y"), textymin, 'SUH 22 senge', fontsize=Fsize,
                 rotation=90, color='gray',
                 horizontalalignment='center', verticalalignment='bottom')

        # --------- LABELS ---------
        plt.xlabel('Dato for måling (kl '+str(measurehour)+')')
        plt.xticks(rotation=30, ha='right')
        plt.ylabel('Belægningsprocent Lungemedicin')

        # --------- RANGES ---------
        plt.ylim(0, plt.gca().get_ylim()[1])
        plt.gca().yaxis.set_major_formatter(PercentFormatter(1*100.))
        # --------- LEGEND ---------
        leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=2, numpoints=1,
                         bbox_to_anchor=(0.47, 1.12), )  # add the legend
        leg.get_frame().set_alpha(0.7)
        # --------------------------

        plt.savefig(plotdir+plotname)
        plt.clf()
        plt.close('all')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'occupancy_movingaverages.pdf'
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

    xvalues = df_results['dates_'+str(measurehour)]
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d-%Y'))
    plt.gca().xaxis.set_major_locator(mdates.YearLocator())
    plt.grid(color='gray', linestyle=':', linewidth=1, alpha=0.5)

    xerr = None
    yerr = None

    plt.grid(linestyle=':', linewidth=lthick/2.)
    # --------- RANGES ---------
    ymin = 0.0
    ymax = 130.0
    dy   = ymax - ymin
    plt.ylim([ymin,ymax])

    # --------- COLORMAP ---------
    cmap = plt.cm.get_cmap('viridis')
    #cmap = plt.cm.get_cmap('plasma')
    cmin = 0
    cmax = 100

    colnorm = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax)
    cmaparr = np.linspace(cmin, cmax, num=50)
    m = plt.cm.ScalarMappable(cmap=cmap)
    m.set_array(cmaparr)
    #cb = plt.colorbar(m)
    #cb.set_label('belægningsprocent')

    # --------- POINT AND CURVES ---------
    colvallist = np.arange(cmin, cmax, cmax/len(measurehours))
    for mm, measurehour in enumerate(measurehours):
        pointcolor = cmap(colnorm(colvallist[mm]))
        plt.errorbar(xvalues, df_results['occupancy_available_movingavg_'+str(measurehour)], xerr=xerr, yerr=yerr,
                     marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                     markerfacecolor=pointcolor, ecolor=pointcolor,
                     markeredgecolor=pointcolor, zorder=20.,
                     label='30 dage glidende gennemsnit (kl '+str(measurehour)+')')

    plt.plot(xvalues, np.zeros(len(xvalues)) + 100, '--', color='black', lw=lthick, zorder=5)
    plt.plot(xvalues, np.zeros(len(xvalues)) + 85, ':', color='black', lw=lthick, zorder=5)

    lineymin = ymin + dy*0.24
    lineymax = ymin + dy*0.40
    textymin = ymin + dy*0.05
    textymax = ymin + dy*0.10
    plt.plot([datetime.datetime.strptime("10-06-2021", "%d-%m-%Y"), datetime.datetime.strptime("10-06-2021", "%d-%m-%Y")],
             [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
    plt.text(datetime.datetime.strptime("10-06-2021", "%d-%m-%Y"), textymin, '10-06-2021', fontsize=Fsize, rotation=90, color='gray',
             horizontalalignment='center', verticalalignment='bottom')

    plt.plot([datetime.datetime.strptime("01-03-2021", "%d-%m-%Y"), datetime.datetime.strptime("01-03-2021", "%d-%m-%Y")],
             [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
    plt.text(datetime.datetime.strptime("01-03-2021", "%d-%m-%Y"), textymin, '01-03-2021', fontsize=Fsize, rotation=90, color='gray',
             horizontalalignment='center', verticalalignment='bottom')

    # --------- LABELS ---------
    plt.xlabel('Dato for måling')
    plt.ylabel('Belægningsprocent Lungemedicin Næstved')

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=2, numpoints=1,
                     bbox_to_anchor=(0.5, 1.12), )  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------

    plt.savefig(plotdir+plotname)
    plt.clf()
    plt.close('all')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'visitationcounts.pdf'
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
    #plt.title('Henvisninger til lungemed. NAE',fontsize=Fsize)

    xvalues = df_results['dates_'+str(measurehour)]
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d-%Y'))
    plt.gca().xaxis.set_major_locator(mdates.YearLocator())

    xerr = None
    yerr = None

    plt.grid(linestyle=':', linewidth=lthick/2.)
    # --------- RANGES ---------
    ymin = -5.0
    ymax = 10.0
    dy   = ymax - ymin
    plt.ylim([ymin,ymax])

    # --------- COLORMAP ---------
    cmap = plt.cm.get_cmap('viridis')
    #cmap = plt.cm.get_cmap('plasma')
    cmin = 0
    cmax = 100

    colnorm = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax)
    cmaparr = np.linspace(cmin, cmax, num=50)
    m = plt.cm.ScalarMappable(cmap=cmap)
    m.set_array(cmaparr)
    #cb = plt.colorbar(m)
    #cb.set_label('belægningsprocent')

    # --------- POINT AND CURVES ---------
    pointcolor = cmap(colnorm(10))
    plt.errorbar(xvalues, df_results['count_vis_aka_23'], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20.,
                 label='Akut afdelingen')

    pointcolor = cmap(colnorm(35))
    plt.errorbar(xvalues, df_results['count_vis_lungNAE_23'], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20.,
                 label='Lungemed NAE (egen afd.)')

    pointcolor = cmap(colnorm(60))
    plt.errorbar(xvalues, df_results['count_vis_lungSLA_23'], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20.,
                 label='Lungemde. SLA')

    pointcolor = cmap(colnorm(85))
    plt.errorbar(xvalues, df_results['count_vis_other_23'], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20.,
                 label='Andre afdelinger')

    lineymin = ymin + dy*0.24
    lineymax = ymin + dy*0.40
    textymin = ymin + dy*0.05
    textymax = ymin + dy*0.10
    plt.plot([datetime.datetime.strptime("10-06-2021", "%d-%m-%Y"), datetime.datetime.strptime("10-06-2021", "%d-%m-%Y")],
             [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
    plt.text(datetime.datetime.strptime("10-06-2021", "%d-%m-%Y"), textymin, '10-06-2021', fontsize=Fsize, rotation=90, color='gray',
             horizontalalignment='center', verticalalignment='bottom')

    plt.plot([datetime.datetime.strptime("01-03-2021", "%d-%m-%Y"), datetime.datetime.strptime("01-03-2021", "%d-%m-%Y")],
             [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
    plt.text(datetime.datetime.strptime("01-03-2021", "%d-%m-%Y"), textymin, '01-03-2021', fontsize=Fsize, rotation=90, color='gray',
             horizontalalignment='center', verticalalignment='bottom')

    # --------- LABELS ---------
    plt.xlabel('Dato for måling')
    plt.ylabel('Henvisninger til Lungemedicin Næstved')

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=2, numpoints=1,
                     bbox_to_anchor=(0.5, 1.12), )  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------

    plt.savefig(plotdir+plotname)
    plt.clf()
    plt.close('all')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = 'visitationcounts_movingsum.pdf'
    if verbose: print(' - Initiating ' + plotname)

    fig = plt.figure(figsize=(9, 6))
    fig.subplots_adjust(wspace=0.1, hspace=0.1, left=0.1, right=0.97, bottom=0.10, top=0.85)
    Fsize = 10
    lthick = 2
    marksize = 4
    plt.rc('text', usetex=False)
    plt.rc('font', family='serif', size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    # plt.title('Henvisninger til lungemed. NAE',fontsize=Fsize)

    xvalues = df_results['dates_' + str(measurehour)]
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d-%Y'))
    plt.gca().xaxis.set_major_locator(mdates.YearLocator())
    plt.grid(color='gray', linestyle=':', linewidth=2)

    xerr = None
    yerr = None

    plt.grid(linestyle=':', linewidth=lthick / 2.)
    # --------- COLORMAP ---------
    cmap = plt.cm.get_cmap('viridis')
    # cmap = plt.cm.get_cmap('plasma')
    cmin = 0
    cmax = 100

    colnorm = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax)
    cmaparr = np.linspace(cmin, cmax, num=50)
    m = plt.cm.ScalarMappable(cmap=cmap)
    m.set_array(cmaparr)
    # cb = plt.colorbar(m)
    # cb.set_label('belægningsprocent')

    # --------- CALC MOVING SUM ---------
    Ndayssum = 30
    df_results['count_vis_aka_23_movingsum'] = df_results['count_vis_aka_23'].rolling(window=Ndayssum).sum()
    df_results['count_vis_lungNAE_23_movingsum'] = df_results['count_vis_lungNAE_23'].rolling(window=Ndayssum).sum()
    df_results['count_vis_lungSLA_23_movingsum'] = df_results['count_vis_lungSLA_23'].rolling(window=Ndayssum).sum()
    df_results['count_vis_other_23_movingsum'] = df_results['count_vis_other_23'].rolling(window=Ndayssum).sum()

    if not concatenateLPR:
        df_SLAupdates['count_vis_aka_23_movingsum']     = df_SLAupdates['count_vis_aka_23'].rolling(window=Ndayssum).sum()
        df_SLAupdates['count_vis_lungNAE_23_movingsum'] = df_SLAupdates['count_vis_lungNAE_23'].rolling(window=Ndayssum).sum()
        df_SLAupdates['count_vis_lungSLA_23_movingsum'] = df_SLAupdates['count_vis_lungSLA_23'].rolling(window=Ndayssum).sum()
        df_SLAupdates['count_vis_other_23_movingsum']   = df_SLAupdates['count_vis_other_23'].rolling(window=Ndayssum).sum()

        df_SUHupdates['count_vis_aka_23_movingsum']     = df_SUHupdates['count_vis_aka_23'].rolling(window=Ndayssum).sum()
        df_SUHupdates['count_vis_lungNAE_23_movingsum'] = df_SUHupdates['count_vis_lungNAE_23'].rolling(window=Ndayssum).sum()
        df_SUHupdates['count_vis_lungSLA_23_movingsum'] = df_SUHupdates['count_vis_lungSLA_23'].rolling(window=Ndayssum).sum()
        df_SUHupdates['count_vis_other_23_movingsum']   = df_SUHupdates['count_vis_other_23'].rolling(window=Ndayssum).sum()

    # --------- RANGES ---------
    ymin = 0.0
    if Ndayssum == 90:
        ymax = 250.0
    else:
        ymax = 100.0
    dy = ymax - ymin
    plt.ylim([ymin, ymax])

    # --------- POINT AND CURVES ---------
    # -------------- BASELINE --------------
    pointcolor = cmap(colnorm(10))
    plt.errorbar(xvalues, df_results['count_vis_aka_23_movingsum'], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20.,
                 label='Akut afdelingen (glidende sum over '+str(Ndayssum)+' dage)')

    pointcolor = cmap(colnorm(35))
    plt.errorbar(xvalues, df_results['count_vis_lungNAE_23_movingsum'], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20.,
                 label='Lungemed NAE (egen afd.; glidende sum over '+str(Ndayssum)+' dage)')

    pointcolor = cmap(colnorm(60))
    plt.errorbar(xvalues, df_results['count_vis_lungSLA_23_movingsum'], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20.,
                 label='Lungemed. SLA (glidende sum over '+str(Ndayssum)+' dage)')

    pointcolor = cmap(colnorm(85))
    plt.errorbar(xvalues, df_results['count_vis_other_23_movingsum'], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20.,
                 label='Andre afdelinger (glidende sum over '+str(Ndayssum)+' dage)')

    #-------------- SLA --------------
    pointcolor = cmap(colnorm(10))
    plt.errorbar(df_SLAupdates['dates_' + str(measurehour)], df_SLAupdates['count_vis_aka_23_movingsum'], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20., linestyle=':',
                 label='SLA LPR (siden 10-03-2022)')

    pointcolor = cmap(colnorm(35))
    plt.errorbar(df_SLAupdates['dates_' + str(measurehour)], df_SLAupdates['count_vis_lungNAE_23_movingsum'], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20., linestyle=':')

    pointcolor = cmap(colnorm(60))
    plt.errorbar(df_SLAupdates['dates_' + str(measurehour)], df_SLAupdates['count_vis_lungSLA_23_movingsum'], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20., linestyle=':')

    pointcolor = cmap(colnorm(85))
    plt.errorbar(df_SLAupdates['dates_' + str(measurehour)], df_SLAupdates['count_vis_other_23_movingsum'], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20., linestyle=':')

    # -------------- SUH --------------
    pointcolor = cmap(colnorm(10))
    plt.errorbar(df_SUHupdates['dates_' + str(measurehour)], df_SUHupdates['count_vis_aka_23_movingsum'], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20., linestyle='--',
                 label='SUH LPR (siden 10-03-2022)')

    pointcolor = cmap(colnorm(35))
    plt.errorbar(df_SUHupdates['dates_' + str(measurehour)], df_SUHupdates['count_vis_lungNAE_23_movingsum'], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20., linestyle='--')

    pointcolor = cmap(colnorm(60))
    plt.errorbar(df_SUHupdates['dates_' + str(measurehour)], df_SUHupdates['count_vis_lungSLA_23_movingsum'], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20., linestyle='--')

    pointcolor = cmap(colnorm(85))
    plt.errorbar(df_SUHupdates['dates_' + str(measurehour)], df_SUHupdates['count_vis_other_23_movingsum'], xerr=xerr, yerr=yerr,
                 marker='.', lw=lthick, markersize=0, alpha=1.0, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20., linestyle='--')

    lineymin = ymin + dy * 0.24
    lineymax = ymin + dy * 0.40
    textymin = ymin + dy * 0.05
    textymax = ymin + dy * 0.10
    # plt.plot(
    #     [datetime.datetime.strptime("10-06-2021", "%d-%m-%Y"), datetime.datetime.strptime("10-06-2021", "%d-%m-%Y")],
    #     [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
    # plt.text(datetime.datetime.strptime("10-06-2021", "%d-%m-%Y"), textymin, '10-06-2021', fontsize=Fsize, rotation=90,
    #          color='gray',
    #          horizontalalignment='center', verticalalignment='bottom')
    #
    # plt.plot(
    #     [datetime.datetime.strptime("01-03-2021", "%d-%m-%Y"), datetime.datetime.strptime("01-03-2021", "%d-%m-%Y")],
    #     [lineymin, lineymax], '-', color='gray', lw=lthick, zorder=5)
    # plt.text(datetime.datetime.strptime("01-03-2021", "%d-%m-%Y"), textymin, '01-03-2021', fontsize=Fsize, rotation=90,
    #          color='gray',
    #          horizontalalignment='center', verticalalignment='bottom')

    # --------- LABELS ---------
    plt.xlabel('Dato for måling')
    plt.ylabel('Henvisninger til Lungemedicin')

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=2, numpoints=1,
                     bbox_to_anchor=(0.48, 1.2), )  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------

    plt.savefig(plotdir + plotname)
    plt.clf()
    plt.close('all')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #Figure plotting number of visitations as a function of bed occupancy; trends between these?
    plotname = 'visitationsVSoccupancy_available.pdf'
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

    xerr = None
    yerr = None

    plt.grid(linestyle=':', linewidth=lthick/2.)
    # --------- RANGES ---------
    ymin = 55.0
    ymax = 135.0
    dy   = ymax - ymin
    plt.ylim([ymin,ymax])

    #plt.xscale('log')
    #plt.yscale('log')

    # --------- COLORMAP ---------
    cmap = plt.cm.get_cmap('viridis')
    #cmap = plt.cm.get_cmap('plasma')
    cmin = 0
    cmax = 100

    colnorm = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax)
    cmaparr = np.linspace(cmin, cmax, num=50)
    m = plt.cm.ScalarMappable(cmap=cmap)
    m.set_array(cmaparr)
    #cb = plt.colorbar(m)
    #cb.set_label('belægningsprocent')

    # --------- POINT AND CURVES ---------
    pointcolor  = cmap(colnorm(10))
    measurehour = 23
    plt.errorbar(df_results['count_vis_aka_23_movingsum'],     df_results['occupancy_actual_movingavg_'+str(measurehour)],
                 xerr=xerr, yerr=yerr,
                 marker='o', lw=0, markersize=marksize, alpha=0.5, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20.,
                 label='Akut afdelingen (glidende sum over '+str(Ndayssum)+' dage)')

    pointcolor = cmap(colnorm(35))
    plt.errorbar(df_results['count_vis_lungNAE_23_movingsum'],     df_results['occupancy_actual_movingavg_'+str(measurehour)],
                 xerr=xerr, yerr=yerr,
                 marker='o', lw=0, markersize=marksize, alpha=0.5, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20.,
                 label='Lungemed NAE (egen afd.; glidende sum over '+str(Ndayssum)+' dage)')

    pointcolor = cmap(colnorm(60))
    plt.errorbar(df_results['count_vis_lungSLA_23_movingsum'],     df_results['occupancy_actual_movingavg_'+str(measurehour)],
                 xerr=xerr, yerr=yerr,
                 marker='o', lw=0, markersize=marksize, alpha=0.5, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=20.,
                 label='Lungemed. SLA (glidende sum over '+str(Ndayssum)+' dage)')

    pointcolor = cmap(colnorm(85))
    plt.errorbar(df_results['count_vis_other_23_movingsum'],     df_results['occupancy_actual_movingavg_'+str(measurehour)],
                 xerr=xerr, yerr=yerr,
                 marker='o', lw=0, markersize=marksize, alpha=0.5, color=pointcolor,
                 markerfacecolor=pointcolor, ecolor=pointcolor,
                 markeredgecolor=pointcolor, zorder=10.,
                 label='Andre afdelinger (glidende sum over '+str(Ndayssum)+' dage)')

    # --------- LABELS ---------
    plt.xlabel('Henvisninger til Lungemedicin Næstved')
    plt.ylabel('Belægningsprocent Lungemedicin Næstved (kl. '+str(measurehour)+')')

    # --------- LEGEND ---------
    leg = plt.legend(fancybox=True, loc='upper center', prop={'size': Fsize / 1.0}, ncol=2, numpoints=1,
                     bbox_to_anchor=(0.5, 1.12), )  # add the legend
    leg.get_frame().set_alpha(0.7)
    # --------------------------

    plt.savefig(plotdir+plotname)
    plt.clf()
    plt.close('all')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Function loading SP data and plotting it against the LPR3 data for comparison.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Predictive time series calculations...
    # https://machinelearningmastery.com/how-to-grid-search-triple-exponential-smoothing-for-time-series-forecasting-in-python/
    # https://machinelearningmastery.com/time-series-forecasting-methods-in-python-cheat-sheet/
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if verbose: print('   Plots saved to \n   ' + plotdir)
#=======================================================================================================================

