"""
#----------------------------
#   NAME
#----------------------------
# emissionlinetests.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# A collection of tests which can be performed to confirm the nature of an
# emission line in a 1D or 2D spectrum.
# Build for testing lines in MOSFIRE data.,
#----------------------------
#   COMMENTS
#----------------------------
# The script contains the following tests:
#     size2D          Testing if the size of the lines is larger than PSF in 2D spec

assymetry(spec1D,linecenter,plot=0): 
size2D(specarr,linecenter,plot=0): 
skyline(lam): 
shadowimg(eps,linecenter):

#
# The wrapper testing4emissionline.py uses these tests to detemine weather a
# line is real or not in reduced MOSFIRE data.

# Tests inspired by Lehnert et al 2010 Nature (so just because the results are positive 
# doesn't mean it's a real line)
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# >>> import emissionlinestests as elt
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-05-03  started by K. B. Schmidt (UCSB)
#----------------------------
"""
#----------------------------
#   MODULES
#----------------------------
import numpy as np         # enable opening with genfromtxt
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scipy.special
import pdb                 # for debugging with pdb.set_trace()
#-------------------------------------------------------------------------------------------------------------
__version__ = 1.0 
__author__ = "K. B. Schmidt (UCSB)"
#-------------------------------------------------------------------------------------------------------------
#                                               TESTS
#-------------------------------------------------------------------------------------------------------------
def assymetry(spec1D,linecenter,plot=0): 
    """
    Determining the size of the line in the epsfitsfile to be able to compare with the psf

    Note: If fit fails with:
          RuntimeError: Optimal parameters not found: Number of calls to function has reached maxfev = 1000.
          try widening the extraction aperture. If the line is offset no peak in the 1D extraction wil
          be present to fit...
    ---- INPUT ----
    Spec1D          The 1D spectral values. Assuming they are on linearlly equidistant wavelength scale 
    linecenter      X-pixel coordiante of center of line

    ---- OUTPUT ----


    ---- EXAMPLE OF USAGE ----

    """
    pixels             = np.arange(len(spec1D))                     # the pixel values corresponding to spectrum

    #spec1D             = gauss(np.arange(len(spec1D)),*[1,linecenter,50]) # test gaussian data

    p0                 = [-1., len(spec1D)/2.0,np.mean(spec1D)]     # initial guess for the fitting coefficients (A, mu and sigma)
    ccont, var_matrix  = curve_fit(parab, pixels, spec1D, p0=p0)    # rough estimate of continuum (or noise profiel if no cont.)
    fitcont            = parab(pixels, *ccont)                      # Get the fitted curve

    p0                 = [1., linecenter, 1.]                       # initial guess for the fitting coefficients (A, mu and sigma)
    cgauss, var_matrix = curve_fit(gauss, pixels, spec1D-fitcont, p0=p0)    # fitting gauss to data
    fitgauss           = gauss(pixels, *cgauss)                     # Get the fitted curve

    if cgauss[2] < 0: print ":: emissionlinetests.assymetry :: WARNING Standard deviation of guassian fit was negative... really?"

    skewent            = range(int(linecenter-3*cgauss[2]),int(linecenter+3*cgauss[2]),1) # pixels +/- 3sigma of line
    #skewent            = range(int(linecenter-3*1.0),int(linecenter+3*1.0),1) # pixels +/- 3sigma of line
    p0                 = np.append(cgauss,0.0)                      # initial guess for the fitting coefficients (A, mu and sigma)
    cskew, var_matrix  = curve_fit(skew, pixels[skewent], spec1D[skewent]-fitcont[skewent], p0=p0)  # fitting skewed gauss to data
    fitskew            = skew(pixels, *cskew)                       # Get the fitted curve

    if plot == 1:
        plt.plot(pixels, spec1D-fitcont, label='Data')
        plt.plot(pixels, fitgauss, label='Fit: Line, Gauss')
        plt.plot(pixels, fitskew, label='Fit: Line, Skew')
        #plt.plot(pixels, fitcont, label='Fit: Continuum')
        leg = plt.legend(fancybox=True, loc='upper right',numpoints=1)
        leg.get_frame().set_alpha(0.6)
        plt.show(block=True)
	

    return cgauss,cskew

#-------------------------------------------------------------------------------------------------------------
def size2D(specarr,linecenter,plot=0): 
    """
    Determines the size of the line in a cutout of the eps fits file.
    Returns the paramters of a least squares fit to a 2D gaussian

    code for 2D fitting inspired by: 
    http://scipy-user.10969.n7.nabble.com/multidimensional-least-squares-fitting-td8456.html

    ---- INPUT ----
    specarr         numpy array with 2D spectrum. E.g. read from the *eps.fits file returned 
                    by the MOSFIRE pipeline.
                    HINT: Use a smaller cutout around line to imporve fit
    linecenter      Pixel coordiantes in spectrum of center of line. Expects [xpix,ypix]
    plot            Set to 1 in order to display plot of results
    ---- OUTPUT ----
    2Dgaussparam    [A,mux,muy,sigmax,sigmay,rho]

    ---- EXAMPLE OF USAGE ----

    """
    shapearr       = specarr.shape
    dataX          = np.arange(shapearr[1])
    dataY          = np.arange(shapearr[0])
    dataXX, dataYY = scipy.meshgrid(dataX,dataY) 
    dataZZ         = specarr
    p0             = [np.max(dataZZ),linecenter[0],linecenter[1],5,5,0.0] # initial guess of [A,mux,muy,sigmax,sigmay,rho]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    test = 0   # Setting up test-case instead of with real data?
    if test == 1:
        dataX          = (np.arange(20)+1)/10.*3
        dataY          = (np.arange(20)+1)/10.*3
        p2D            = [1.0,2.1111111,2.5,0.75,0.5,0.83]
        dataXX, dataYY = scipy.meshgrid(dataX,dataY)     
        dataZZ         = gauss_2Darr(dataXX,dataYY,*p2D)
        p0             = [0.5,2.0,2.2,0.3,0.3,0.0]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    p_res, p_cov = scipy.optimize.leastsq(gauss_2Darr_residual, x0=p0, args=(dataZZ, dataXX, dataYY)) 
    #print p_res

    if plot == 1:
        import pylab 
        from mpl_toolkits.mplot3d import Axes3D
        fig = pylab.figure() 
        ax = Axes3D(fig) 
        ax.plot_wireframe(dataXX, dataYY, dataZZ,label='DATA') 

        ZZfit = gauss_2Darr(dataXX,dataYY,*p_res)
        ax.plot_wireframe(dataXX, dataYY, ZZfit,color='r',label='MODEL') 

        ax.set_xlabel('x') 
        ax.set_ylabel('y')
        ax.set_zlabel('z') 
        ax.legend(fancybox=True, loc='upper left')
        pylab.show(block=True) 

    return p_res, p_cov
#-------------------------------------------------------------------------------------------------------------
def skyline(lam): 
    """
    Checking if there are prominent sky lines at location of line
    ---- INPUT ----
    lam  = range of wavelengths in Angstrom to check for skylines
           given on the form [min(lam),max(lam)]

    ---- OUTPUT ----
    linematch:  numpy array with lines in lam-region.
                linematch[:,0] = wavelentgh of OH line
                linematch[:,1] = strength of OH line
    """
    OHfile  = '/Users/kasperborelloschmidt/work/observing/OHairglowDAT/rousselot2000.dat.txt'
    OHdat   = np.genfromtxt(OHfile,dtype=None,comments='#')
    lineent = np.where((OHdat[:,0] > lam[0]) & (OHdat[:,0] < lam[1]))
    Nmatch  = len(lineent[0])

    if Nmatch > 0:
        linematch = OHdat[lineent]
    else:
        linematch = np.array([[-99.,-99]])

    return linematch
#-------------------------------------------------------------------------------------------------------------
    # EL offset in slit?

#-------------------------------------------------------------------------------------------------------------
    # EL = Noise peak?

#-------------------------------------------------------------------------------------------------------------
def shadowimg(eps,linecenter):
    """
    Looking for the expected two shadow images in the eps data.
    Use picel postion of line, offsets 2.5'' above an below, iverts image, and fit a 2D gaussian.
    if amplitude is a ~0.5 of the actual line a detection of a shadow is declared.
    """



    shadowpos = [1,1,1,1]
    return shadowpos
#-------------------------------------------------------------------------------------------------------------
#                                               UTILITIES
#-------------------------------------------------------------------------------------------------------------
def gauss(x, *p):
    """
    Gaussian/normal distribution
    """
    A, mu, sigma = p
    t = (x-mu)/sigma
    return A*np.exp(-t**2/2.)
#-------------------------------------------------------------------------------------------------------------
def gauss_2D(xyvec, *p):
    """
    --- SEE gauss_2D_arr INSTEAD ---

    Gaussian/normal distribution in 2D
    The parmeters are:
       [A,mux,muy,sigmax,sigmay]
    --- EXAMPLE ---
    x   = [1,2,3]
    y   = [1,2,3]
    p   = [1.0,2.0,2.5,0.25,0.5]
    G2D = gauss_2D(np.array([x,y]), *p)
    """
    A, mux, muy, sigmax, sigmay = p
    x         = xyvec[:,0]
    y         = xyvec[:,1]
    Nx, Ny    = len(x), len(y)
    rho       = 0.0                               # the correlation between x and y where xsigma and ysigma are > 0
    amplitude = A/np.sqrt(1-rho**2.)
    gauss2D   = np.zeros((Nx,Ny),dtype=np.float32)

    for ii in range(Nx):
        for jj in range(Ny):
            tx  = (x[ii]-mux)/sigmax
            ty  = (y[jj]-muy)/sigmay
            gauss2D[ii,jj] = amplitude*np.exp( -1 * (tx**2 + ty**2 - 2*rho*tx*ty) / (2.*(1-rho**2.0)) )

    # could also do: gauss2D = np.multiply(gauss(x, *px).reshape(shapearr[0],1),gauss(y, *py).reshape(1,shapearr[1]))

    return gauss2D
#-------------------------------------------------------------------------------------------------------------
def gauss_2Darr(xx,yy, *p):
    """
    Gaussian/normal distribution in 2D

    The parmeters are:
       [A,mux,muy,sigmax,sigmay,rho]
    where rho is the correlation between x and y where xsigma and ysigma are > 0

    --- EXAMPLE ---
    x   = np.array([1,2,3])
    y   = np.array([1,2,3])
    xx,yy = np.meshgrid(x,y)
    p   = [1.0,2.0,2.5,0.25,0.5]
    G2D = gauss_2Darr(xx,yy, *p)
    """
    A, mux, muy, sigmax, sigmay, rho = p
    amplitude = A/np.sqrt(1-rho**2.)
    tx  = (xx-mux) / sigmax
    ty  = (yy-muy)/ sigmay

    gauss2D = amplitude*np.exp( -1 * (tx*tx + ty*ty - 2*rho*tx*ty) / (2.*(1-rho**2.0)) )

    return gauss2D
#------------------------------------------------------------------------------ 
def gauss_2Darr_residual(p, G2Dinput, xx, yy): 
    """ 
    Residual from Gauss2D function
    """ 
    G2D = gauss_2Darr(xx,yy,*p)
    return np.ravel(G2D - G2Dinput)
#-------------------------------------------------------------------------------------------------------------
def skew(x, *p):
    """
    A skewed normal distribution where alpha=0 is not skewed (i.e. = gauss)
    """
    A, mu, sigma, alpha = p

    t = (x-mu)/sigma

    normpdf = A * np.exp(-t**2/2)
    normcdf = (1 + scipy.special.erf(alpha*t/np.sqrt(2))) / 2.

    return 2.0 / sigma * normpdf * normcdf
#-------------------------------------------------------------------------------------------------------------
def parab(x, *p):
    A, B, C = p
    return A*x**2 + B*x + C

#-------------------------------------------------------------------------------------------------------------
#                                                  END
#-------------------------------------------------------------------------------------------------------------


