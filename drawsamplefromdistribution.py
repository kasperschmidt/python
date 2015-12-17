#=============================================================================================================
#=============================================================================================================
def drawsamplefromdistribution(points1D,bins1D,points2D=None,bins2D=None,verbose=False,size=1,stop=False,plot=False):
    """
    ----------------------------
       NAME
    ----------------------------
     drawsamplefromdistribution.py
    ----------------------------
       PURPOSE/DESCRIPTION
    ----------------------------
     Using the inverted cumulative distribution function (ICDF) method to draw samples
     from a provided distribution of points in 1D or 2D
    ----------------------------
       COMMENTS
    ----------------------------
     NOTE... you might want to use the generalized version in drawsamplefromdistribution_multiD.py
    ----------------------------
       INPUTS:
    ----------------------------
     points1D         : points making up the distribution to draw values from (x dimension).
     bins1D           : number of bins to use for histogram
    ----------------------------
       OPTIONAL INPUTS:
    ----------------------------
     points2D         : points making up the distribution to draw values from in second dimension (y dimension).
     bins2D           : number of bins to use for histogram
     size             : size of sample to draw (DEFAULT is 1) 
     verbose          : set True to get info/messages printed to the screen
     stop             : if True stop at end before returning draws 
     plot             : if True distribution and draws are plotted.
    ----------------------------
       OUTPUTS:
    ----------------------------
     Numpy array containing draws from distribution
    ----------------------------
       EXAMPLES/USAGE
    ----------------------------
     >>> import drawsamplefromdistribution as dfd
     >>> draws = dfd.drawsamplefromdistribution(samplex,binsx,points2D=sampley,bins2D=biny,verbose=True,size=10,plot=True)
    ----------------------------
       BUGS
    ----------------------------
    
    ----------------------------
       REVISION HISTORY
    ----------------------------
     2013-08-08  started by K. B. Schmidt (UCSB)
    ----------------------------
    """
    #-------------------------------------------------------------------------------------------------------------
    import numpy as np         # enable opening with genfromtxt
    import pdb                 # for debugging with pdb.set_trace()
    #-------------------------------------------------------------------------------------------------------------
    # Drawing samples from 1D distribution
    Ndraws    = size
    hist1D    = np.histogram(points1D,bins=bins1D)
    binedge1D = hist1D[1]
    drawval1D = drawfrom1D(points1D,bins1D,Ndraws)
    #-------------------------------------------------------------------------------------------------------------
    # Distinguishing between 1D and 2D case.
    if points2D == None:
        drawval = drawval1D # defining output
    else: # the 2D case
        drawval2D = np.zeros((Ndraws,4)) # array to contain x and y pairs of draws 
        drawval2D[:,0] = drawval1D[:,0]
        drawval2D[:,1] = drawval1D[:,1]

        for jj in xrange(Ndraws):  
            edgediff  = binedge1D - drawval2D[jj,0]
            bent      = np.where(np.abs(edgediff) == np.min(np.abs(edgediff)))[0][0]    # 1D bin edge closest to draw
            if edgediff[bent] >= 0: # right edge of bin to slice
                ent_slice = np.where((points1D > binedge1D[bent-1]) & (points1D < binedge1D[bent]))
            elif edgediff[bent] < 0:# left edge of bin to slice
                ent_slice = np.where((points1D > binedge1D[bent]) & (points1D < binedge1D[bent+1]))
                
            points2D_slice  = points2D[ent_slice] # Slice in 2D distribution at draw on 1D axis
            if len(points2D_slice) != 0:
                drawval2D_slice = drawfrom1D(points2D_slice,bins2D,1)
                drawval2D[jj,2] = drawval2D_slice[:,0]
                drawval2D[jj,3] = drawval2D_slice[:,1]
            else: # if slice is empty set the draw to NaN.
                if verbose: print '   No values in slice for draw ',jj,' setting to NaN'
                drawval2D[jj,2] = np.nan
                drawval2D[jj,3] = np.nan


            #if drawval2D[jj,2] > 8.0: pdb.set_trace()
            #plt.plot(bincen2D,np.cumsum(np.histogram(points2D_slice,bins=binedge2D)[0]/float(len(points2D_slice))))
        drawval = drawval2D #defining output
    #-------------------------------------------------------------------------------------------------------------
    if (plot == True) & (points1D != None):
        import matplotlib.pyplot as plt
        plotname = 'drawsamplefromdistribution_draw1D.pdf'
        if verbose: print '   Creating plot ',plotname
        plt.clf()

        hist1D     = np.histogram(points1D,bins=bins1D)
        binsize1D  = (np.max(points1D)-np.min(points1D))/bins1D
        binedge1D  = hist1D[1]
        bincen1D   = binedge1D[0:-1]+binsize1D/2.

        dist       = hist1D[0]
        distcum    = np.cumsum(dist)                # cumulative distribution function
        distcumN   = distcum / float(len(points1D)) # normalized CDF

        plt.plot(bincen1D,distcumN,'r-',lw=3,label='Cumulative distribution')
        plt.plot(drawval[:,0],drawval[:,1],'ko',label='The '+str(Ndraws)+' draws')

        plt.xlabel(r'Points 1D')
        plt.ylabel(r'CDF')

        leg = plt.legend(fancybox=True, loc='upper left',numpoints=1,prop={'size':12})
        leg.get_frame().set_alpha(0.6)

        plt.savefig(plotname)    
    #-------------------------------------------------------------------------------------------------------------
    if (plot == True) & (points2D != None):
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm

        plotname = 'drawsamplefromdistribution_draw2D.pdf'
        if verbose: print '   Creating plot ',plotname
        plt.clf()
        
        heatmap, xedges, yedges = np.histogram2d(points2D,points1D,bins=(bins2D,bins1D))
        extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]

        plt.imshow(heatmap, extent=extent, interpolation='nearest',norm=LogNorm(vmin=1, vmax=1000))

        goodent = ~np.isnan(drawval[:,2])
        
        plt.plot(drawval[goodent,0],drawval[goodent,2],'k.',label=r'The '+str(Ndraws)+' draws')

        xmin, xmax, ymin, ymax = -10, 10, -10, 10
        plt.xlim(xmin,xmax)
        plt.ylim(ymin,ymax)
        
        plt.xlabel(r'Points 1D')
        plt.ylabel(r'Points 2D')

        cbar = plt.colorbar()
        cbar.set_label('2D Histogram ('+str(bins1D)+'x'+str(bins2D)+') of Points 1D and Points 2D')

        leg = plt.legend(fancybox=True, loc='upper left',numpoints=1,prop={'size':12})
        leg.get_frame().set_alpha(0.6)

        plt.savefig(plotname)    
    #-------------------------------------------------------------------------------------------------------------
    if stop: pdb.set_trace()
    return drawval
#=============================================================================================================
#=============================================================================================================
def drawfrom1D(pts,bins,Ndraws):
    '''
    Drawing from 1D distribution using the inverted cumulative distribution function (ICDF) 
    method to draw samples

    Expects the binedges (so multiple calls have same bins) which can be obtained from
    hist     = np.histogram(pts,bins=bins)
    binedge  = hist[1]

    Returns numpy array (Ndraws,2) containing the 
        Drawn values       drawval[:,0]
        CDF value of draw  drawval[:,1]
    '''
    import numpy as np
    import pdb
    
    hist     = np.histogram(pts,bins=bins)
    binedge  = hist[1]
    binsize  = (np.max(binedge)-np.min(binedge))/bins
    bincen   = binedge[0:-1]+binsize/2.

    dist       = hist[0]
    distcum    = np.cumsum(dist)                # cumulative distribution function
    distcumN   = distcum / float(len(pts)) # normalized CDF
        
    drawsU     = np.random.uniform(low=0.0, high=1.0, size=Ndraws) # uniformly drawn values on y axis of CDF
    drawval  = np.zeros((Ndraws,2)) # array to contain x and y values of draws

    for dd in xrange(Ndraws):
        absdiff = np.abs(distcumN - drawsU[dd])
        minent  = np.where(absdiff == np.min(absdiff))[0]
        minent  = minent[int(len(minent)/2.)] # taking the 'center' value in case of multiple entries
        
        drawval[dd,0] = (minent+1)*binsize+np.min(binedge)  # x values
        drawval[dd,1] = drawsU[dd] # distcumN[minent]                   # y values
    return drawval
#=============================================================================================================
#=============================================================================================================
