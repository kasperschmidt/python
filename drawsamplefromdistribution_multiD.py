#=============================================================================================================
#=============================================================================================================
def drawsample(dist,Nbins,verbose=False,size=1,Ncut=5,stop=False,plot=False):
    """
    ----------------------------
       NAME
    ----------------------------
     drawsample.py
    ----------------------------
       PURPOSE/DESCRIPTION
    ----------------------------
     Using the inverted cumulative distribution function (ICDF) method to draw a sample
     from a provided distribution multi dimensional distribution of points.
    ----------------------------
       COMMENTS
    ----------------------------
     This is the generalized version of drawsamplefromdistribution.py
    ----------------------------
       INPUTS:
    ----------------------------
     dist             : Array containing the values of each N dimesions of the distribution.
                        Expects a numpy array of size (N,M) where N is the number of dimensions
                        and M is the number of points in the ND distribution
     Nbins            : The number of bins to use in each dimension. I.e. the Nth dimension will
                        be turned into a histogram with the bins 
                        np.linspace(min(points[N,:]),max(points[N,:]),Nbins+1)
    ----------------------------
       OPTIONAL INPUTS:
    ----------------------------
     Ncut             : minimum required number of objects in distribution slice to make a draw
                        Only effective for Ndimensions > 1 as there are no conditional distributions
                        in 1D
     size             : size of sample to draw (DEFAULT is 1) 
     verbose          : set True to get info/messages printed to the screen
    ----------------------------
       OUTPUTS:
    ----------------------------
     Numpy array containing draws from distribution
    ----------------------------
       EXAMPLES/USAGE
    ----------------------------
import numpy as np
import drawsamplefromdistribution_multiD as dfd
Hval, Jval  = np.random.randn(300000), np.random.randn(300000) 
dist  = np.asarray([Jval,Hval])
Nbins = 300
draws = dfd.drawsample(dist,Nbins,verbose=True,size=1000)
Ndim  = dist.shape[0]
d1, d2  = 0, 1
dfd.plot1D(dist[d1,:],draws[:,d1*Ndim+0:d1*Ndim+2],Nbins,save='drawsamplefromdistributionMD_draw1D.pdf')
dfd.plot2D(dist[d1,:],dist[d2,:],draws[:,d1*Ndim+0:d1*Ndim+2],draws[:,d2*Ndim+0:d2*Ndim+2],Nbins,save='drawsamplefromdistributionMD_draw2D.pdf')
    ----------------------------
       BUGS
    ----------------------------
    Note really a bug but note that deviations occur for small number of bins and small samples. 
    ----------------------------
       REVISION HISTORY
    ----------------------------
     2013-08-28  started by K. B. Schmidt (UCSB)
    ----------------------------
    """
    #-------------------------------------------------------------------------------------------------------------
    import numpy as np         # enable opening with genfromtxt
    import pdb                 # for debugging with pdb.set_trace()
    #-------------------------------------------------------------------------------------------------------------
    if   len(dist.shape) == 1:
        Ndim   = 1
        dist1D = dist
    elif len(dist.shape) == 2:
        Ndim   = dist.shape[0]
        dist1D = dist[0,:]
    else:
        print ':: drawsample :: ERROR - shape of distributuon array has neither length 1 or 2 --> ABORTING'
    Ndraws    = size
    drawval   = np.zeros((Ndraws,2*Ndim)) # array to contain x and y pairs of draws 
    #-------------------------------------------------------------------------------------------------------------
    # draw for the first dimension
    draws          = drawfrom1D(dist1D,Nbins,Ndraws)
    drawval[:,0:2] = draws                                                              # fill output array
    #-------------------------------------------------------------------------------------------------------------
    # draw for remaining dimensions
    for jj in xrange(Ndraws):                                                           # looping over number of draws
        binedgeprev    = np.linspace(np.min(dist1D),np.max(dist1D),Nbins+1)             # bin edges of first dimesion
        dimdist        = dist                                                           # resetting dimdist
        #print 'resetting : ',dimdist.shape
        for ii in np.arange(Ndim-1)+1:                                                  # looping over dimensions
            if np.isnan(drawval[jj,(ii-1)*2]):     # if previous slice was empty set this draw to NaN as well
                drawval[jj,ii*2+0] = np.nan
                drawval[jj,ii*2+1] = np.nan
                continue
            edgediff  = binedgeprev - drawval[jj,(ii-1)*2]                              # locating position of draw
            bent      = np.where(np.abs(edgediff) == np.min(np.abs(edgediff)))[0][0]    # bin edge closest to draw
            
            if edgediff[bent] >= 0: # right edge of bin to slice
                ent_slice = np.where((dimdist[ii-1,:] > binedgeprev[bent-1]) & (dimdist[ii-1,:] < binedgeprev[bent]))
            elif edgediff[bent] < 0:# left edge of bin to slice
                try:
                    ent_slice = np.where((dimdist[ii-1,:] > binedgeprev[bent]) & (dimdist[ii-1,:] < binedgeprev[bent+1]))
                except:
                    pdb.set_trace()
                    
            dist_slice  = dimdist[ii,ent_slice].transpose()  # Slice in current distribution for previous draw
            if len(dist_slice) >= Ncut:                      # checking that there are at least 5 points left in slice
                draws_slice = drawfrom1D(dist_slice,Nbins,1)
                drawval[jj,ii*2+0] = draws_slice[:,0]
                drawval[jj,ii*2+1] = draws_slice[:,1]
            else:                                         # if slice is empty set the draw to NaN.
                if verbose: print '   < ',Ncut,' values in slice of dimension ',ii+1,' for draw ',jj,' --> setting to NaN'
                drawval[jj,ii*2+0] = np.nan
                drawval[jj,ii*2+1] = np.nan

            try:
                binedgeprev = np.linspace(np.min(dimdist[ii,:]),np.max(dimdist[ii,:]),Nbins+1)  # updating binedgeprev
                dimdist     = dimdist[:,ent_slice[0]]                                           # updating dimdist 
                if np.max(binedgeprev) < drawval[jj,ii*2+0]: pdb.set_trace()                    # error checking
            except:
                if ~np.isnan(drawval[jj,ii*2+0]): pdb.set_trace()                               # error checking
                continue # if dimdist is empty move on (the rest of the draws will be nans
        #if jj == 10: pdb.set_trace()
    #-------------------------------------------------------------------------------------------------------------    
    return drawval
#=============================================================================================================
#=============================================================================================================
def drawfrom1D(pts,bins,Ndraws):
    """
    Drawing from 1D distribution using the inverted cumulative distribution function (ICDF) 
    method to draw samples

    Expects the binedges (so multiple calls have same bins) which can be obtained from
    hist     = np.histogram(pts,bins=bins)
    binedge  = hist[1]

    Returns numpy array (Ndraws,2) containing the 
        Drawn values       drawval[:,0]
        CDF value of draw  drawval[:,1]
    """
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
        
        drawval[dd,0] = (minent+1)*binsize+np.min(binedge)-binsize/2.   # x values
        drawval[dd,1] = drawsU[dd] # distcumN[minent]                   # y values
    return drawval
#=============================================================================================================
#=============================================================================================================
def plot1D(dim1,draw1,Nbins,xlab='Points 1D',save=False):
    """
    Plotting a 1 dimentiosnal draw (commparing with input)
    To save the figure instead of showing it give desired output name to save keyword 
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import pdb
    plt.clf()
    goodent    = ~np.isnan(draw1[:,0])
    Ndraws     = draw1[goodent,:].shape[0]

    hist1D     = np.histogram(dim1,bins=Nbins)
    binsize1D  = (np.max(dim1)-np.min(dim1))/Nbins
    binedge1D  = hist1D[1]
    bincen1D   = binedge1D[0:-1]+binsize1D/2.

    dist       = hist1D[0]
    distcum    = np.cumsum(dist)                # cumulative distribution function
    distcumN   = distcum / float(len(dim1))     # normalized CDF

    plt.plot(bincen1D,distcumN,'r-',lw=3,label='Cumulative distribution')
    plt.plot(draw1[goodent,0],draw1[goodent,1],'ko',label='The '+str(Ndraws)+' draws (!= NaN)')

    plt.xlabel(xlab)
    plt.ylabel(r'CDF')

    leg = plt.legend(fancybox=True, loc='upper left',numpoints=1,prop={'size':12})
    leg.get_frame().set_alpha(0.6)

    if save == False: 
        plt.show()
    else:
        plt.savefig(save)    
#=============================================================================================================
#=============================================================================================================
def plot2D(dim1,dim2,draw1,draw2,Nbins,xlab='Points 1D',ylab='Points 2D',save=False):
    """
    Plotting a 2 dimentiosnal draw (commparing with input)
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import numpy as np
    import pdb
    plt.clf()
    goodent = ( ~np.isnan(draw1[:,0]) ) & ( ~np.isnan(draw2[:,0]) )
    Ndraws     = draw1[goodent,:].shape[0]
        
    heatmap, xedges, yedges = np.histogram2d(dim2,dim1,bins=(Nbins,Nbins))
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]

    plt.imshow(heatmap, extent=extent, interpolation='nearest',norm=LogNorm(vmin=1, vmax=1000))

    plt.plot(draw1[goodent,0],draw2[goodent,0],'k.',label=r'The '+str(Ndraws)+' draws (!= NaN)')

    xmin, xmax, ymin, ymax = -10, 10, -10, 10
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
        
    plt.xlabel(xlab)
    plt.ylabel(ylab)

    cbar = plt.colorbar()
    cbar.set_label('2D Histogram ('+str(Nbins)+'x'+str(Nbins)+') of '+xlab+' and '+ylab)

    leg = plt.legend(fancybox=True, loc='upper left',numpoints=1,prop={'size':12})
    leg.get_frame().set_alpha(0.6)

    if save == False: 
        plt.show()
    else:
        plt.savefig(save)    
#=============================================================================================================
#=============================================================================================================
