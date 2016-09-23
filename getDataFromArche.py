# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import commands
import glob
import pyfits
import numpy as np
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def download_data(archeuser,field='cosmos',pointing=10,collection='QtClassify',outputdir='fielddir',
                  port='2222',acsimg='814w',acsimgvs='1.0',lsdcatvs='1.0',SNstr='',download=True,clobber=False,verbose=True):
    """

    Downloading data for a given MUSE-Wide pointing via SCP

    --- INPUT ---
    field            the field to download data for. 'cosmos' or 'cdfs'
    pointing         the MUSE-Wide pointing number to grab data for.
    type             Specify the collection of data to download. Choose between
                      'QtCLassify'    Data products needed to run QtClassify
                      'all'           All products found in field+pointing directory
    outputdir        The output directory to save the downloaded files to.
                     By default ('fielddir') the files will be stored in the ./candels-*field*-*pointing*/
                     directory assuming it exisits.
    port             Port to use when connecting to arche
    acsimg           Specify ACS image cutout to download           (used for collection='QtClassify')
    acsimgvs         Version of acs image cutours to download       (used for collection='QtClassify')
    lsdcatvs         specify LSDCat version to download files for   (used for collection='QtClassify')
    SNstr            if LSDCat catalog (cat*) was appended SN string provide it here, .e.g, SNstr='_sn5.0'
    download         If true the data will be downloaded. Otherwise the commands used
                     to download the data will just be printed to the screen.
    clobber          Overwrtite files in output directory if they already exist?
    verbose          Toggle verbosity

    --- EXAMPLE OF USE ---
    import getDataFromArche as gd
    user     = 'jondoe'
    filelist = gd.download_data(user,field='cosmos',pointing='10',outputdir='temp/',collection='all')

    filelist = gd.download_data(user,field='cdfs',pointing='14',collection='QtClassify',SNstr='_sn5.0',lsdcatvs='2.0',acsimgvs='1.0',clobber=False)

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dirname  = 'candels-'+field+'-'+str(pointing)
    if outputdir == 'fielddir':
        outputdir = './'+dirname+'/'
    if verbose: print ' - Will download files from '+dirname+\
                      ' collecting files for the collection='+"'"+collection+"'"

    if verbose: print ' - Putthing together scp command and setting up file lists'
    basedcmd = 'scp -P '+str(port)+' '+archeuser+'@arche.aip.de:/store/data/musewide/'+dirname+'/'

    if collection == 'all':
        lscmd    = 'ssh -p '+str(port)+' '+archeuser+'@arche.aip.de  ls /store/data/musewide/'+dirname+'/*.*'
        lsout    = commands.getoutput(lscmd)

        filesALL = lsout.split('\n')
        filelist = ['*.*']
    elif collection== 'QtClassify':
        filelist = []
        filelist.append('median_filtered_DATACUBE_'+dirname+'_v1.0.fits')
        filelist.append('s2n_opt_v250_'+dirname+'_v'+lsdcatvs+'.fits')
        filelist.append('cat_opt_v250_'+dirname+'_v'+lsdcatvs+SNstr+'.fits')
        filelist.append('acs_'+acsimg+'_'+dirname+'_cut_v'+acsimgvs+'.fits')
    else:
        if verbose: print " - WARNING didn't recognize the collection="+collection+" so returning empty list "
        return []
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Spawning the following commands to the shell:'
    filecounter = 0
    skipcounter = 0
    for archefile in filelist:
        scpcmd = basedcmd+archefile+' '+outputdir

        if verbose: print '   '+scpcmd
        if (clobber == False) & ( len(glob.glob(outputdir+'/'+archefile)) != 0):
            if verbose: print '   file already exists in output directory and clobber=False so moving on'
            skipcounter = skipcounter + 1
        else:
            if download:
                scpout = commands.getoutput(scpcmd)

                if scpout == '':
                    filecounter = filecounter + 1
                else:
                    print scpout

    if collection == 'all':
        filecounter = len(filesALL)
        filelist    = filesALL
    if download:
        if verbose: print ' - Succesfullt downloaded '+str(filecounter)+' / '+str(len(filelist)-skipcounter)+\
                          ' (skipping '+str(skipcounter)+') files from arche '
    else:
        if verbose: print ' - Download=False so no files downloaded from arche'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if collection == 'QtClassify':
        if verbose:
            datacube = filelist[0]
            LSDCatSN = filelist[1]
            LSDCat   = filelist[2]
            HSTimg   = filelist[3]

            print '\n - To run QtClassify move to outputdir ('+outputdir+') and execute (in your shell):'
            print """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    datapath='%s'
    datacube=$datapath'%s'
    LSDCatSN=$datapath'%s'
    LSDCat=$datapath'%s'
    HSTimg=$datapath'%s'
    output=$datapath'%s_QtClassify_output_RENAME_.fits'

    qtclassify -id $datacube -isn $LSDCatSN -c $LSDCat -o $output -F 0 -N 2 -hst $HSTimg --replaceCubeNaNs False --column_X X_PEAK_SN --column_Y Y_PEAK_SN --column_Z Z_PEAK_SN --column_RA RA_PEAK_SN --column_DEC DEC_PEAK_SN --column_LAM LAMBDA_PEAK_SN

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            """ % (outputdir,datacube,LSDCatSN,LSDCat,HSTimg,dirname)
            print '   (here "qtclassify" is an alias for "python' \
                  ' /Local/Path/To/qtclassify/line_classification_GUI_pyqtgraph.py")'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return filelist
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def summarize_QtClassifyOutput(qtclassifyoutputfile,idcol='ID',wavecol='LAMBDA_PEAK_SN',SNcol='DETSN_MAX',verbose=True):
    """

    Generate an ascii file with a condensed summary of a QtClassfication fits output

    --- INPUT ---
    qtclassifyoutputfile        output file from QtCLassify classification to summarize/condense
    idcol                       column name of ID column in fits file
    wavecol                     column name of wavelength column in fits file
    SNcol                       column name of signal to noise of detected lines
    verbose                     toggle verbosity

    --- OUPTUT ---
    summaryfile                 The file name of the ascii file containing the summary generated

    --- EXAMPLE OF USE ---
    import getDataFromArche as gd
    path        = '/Users/kschmidt/work/MUSE/QtClassify/candels-cdfs-18/'
    QtCfile     = path+'candels-cdfs-18_training_QtClassify_output_kschmidt160922.fits'
    summaryfile = gd.summarize_QtClassifyOutput(QtCfile)

    """
    classdat = pyfits.open(qtclassifyoutputfile)[1].data

    outputfile = qtclassifyoutputfile.replace('.fits','_summary.txt')

    fout = open(outputfile,'w')
    fout.write('# Condensed Summary of '+qtclassifyoutputfile+'\n')
    fout.write('# Showing results for strongest line. \n')
    fout.write('# identification: Cres = Continuum residual; Rev = revisit \n')
    fout.write("# id SN_strongestline wavelength Nlines identification redshift quality confidence association #LW# all line wavelengths #C# comments \n")

    IDs  = np.sort(np.unique(classdat[idcol]))
    Nobj = len(IDs)
    if verbose: print ' - Found '+str(Nobj)+' objects in:\n   '+qtclassifyoutputfile+'\n   to summarize classifications for'

    for objid in IDs:
        objent   = np.where(classdat[idcol] == objid)[0]
        Nlines   = len(objent)
        objdat   = classdat[objent]
        if Nlines == 1:
            maxSNent = 0
            short    = objdat['Short'][maxSNent]
            comments = objdat['Comment'][maxSNent]
            objQ     = objdat['Quality'][maxSNent]
            objC     = objdat['Confidence'][maxSNent]
            assoc    = objdat['Association'][maxSNent]
        else:
            maxSNent = np.where(objdat[SNcol] == np.max(objdat[SNcol]))[0]
            short    = objdat['Short'][maxSNent][0]
            comments = objdat['Comment'][maxSNent][0]
            objQ     = objdat['Quality'][maxSNent][0]
            objC     = objdat['Confidence'][maxSNent][0]
            assoc    = objdat['Association'][maxSNent][0]


        lineSN   = objdat[SNcol][maxSNent]
        linewave = objdat[wavecol][maxSNent]

        objz     = objdat['Redshift'][maxSNent]


        linelams = ','.join([str("%.f" % lw) for lw in objdat[wavecol]])

        objstring = str("%.4d" % objid)+'  '+str("%7.2f" % lineSN)+'  '+str("%.f" % linewave)+'  '+str("%.4d" % Nlines)+\
                    ' '+str("%6s" % short)+'   '+str("%.5f" % objz)+'   '+str(objQ)+' '+' '+str(objC)+'   '+\
                    str("%6s" % assoc)+' #O# '+str(linelams)+' #C# '+str(comments)+' \n'

        fout.write(objstring)
    fout.close()
    if verbose: print ' - wrote summary to:\n   '+outputfile

    #import pdb; pdb.set_trace()
    return outputfile
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =