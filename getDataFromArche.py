# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import commands
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def download_data(archeuser,field='cosmos',pointing=10,collection='QtClassify',outputdir='./',
                  port='2222',acsimg='606w',lsdcatvs='1.0',download=True,verbose=True):
    """

    Downloading data for a given MUSE-Wide pointing via SCP

    --- INPUT ---
    field            the field to download data for. 'cosmos' or 'cdfs'
    pointing         the MUSE-Wide pointing number to grab data for.
    type             Specify the collection of data to download. Choose between
                      'QtCLassify'    Data products needed to run QtClassify
                      'all'           All products found in field+pointing directory
    outputdir        The output directory to save the downloaded files to
    port             Port to use when connecting to arche
    acsimg           Specify ACS image cutout to download           (used for collection='QtClassify')
    lsdcatvs         specify LSDCat version to download files for   (used for collection='QtClassify')
    download         If true the data will be downloaded. Otherwise the commands used
                     to download the data will just be printed to the screen.
    verbose          Toggle verbosity

    --- EXAMPLE OF USE ---
    import getDataFromArche as gd
    user     = 'jondoe'
    filelist = gd.download_data(user,field='cosmos',pointing='10',outputdir='temp/',collection='all')

    filelist = gd.download_data(user,field='cdfs',pointing='04',outputdir='temp/',collection='QtClassify')

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dirname  = 'candels-'+field+'-'+str(pointing)
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
        filelist.append('median_filtered_DATACUBE_'+dirname+'_v'+lsdcatvs+'.fits')
        filelist.append('s2n_opt_v250_'+dirname+'_v1.0.fits')
        filelist.append('cat_opt_v250_'+dirname+'_v1.0.fits')
        filelist.append('acs_'+acsimg+'_'+dirname+'_cut.fits')
    else:
        if verbose: print " - WARNING didn't recognize the collection="+collection+" so returning empty list "
        return []
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Spawning the following commands to the shell:'
    filecounter = 0
    for archefile in filelist:
        scpcmd = basedcmd+archefile+' '+outputdir

        if verbose: print '   '+scpcmd
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
        if verbose: print ' - Succesfullt downloaded '+str(filecounter)+' / '+str(len(filelist))+' files from arche '
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
    output=$datapath'QtClassify_output_RENAME_.fits'

    qtclassify -id $datacube -isn $LSDCatSN -c $LSDCat -o $output -F 0 -N 2 -hst $HSTimg --replaceCubeNaNs False
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            """ % (outputdir,datacube,LSDCatSN,LSDCat,HSTimg)
            print '   (here "qtclassify" is an alias for "python' \
                  ' /Local/Path/To/qtclassify/line_classification_GUI_pyqtgraph.py")'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return filelist
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =