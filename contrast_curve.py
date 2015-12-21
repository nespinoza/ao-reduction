import numpy as np
import pyfits
import Utils

def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
        xo = float(xo)
        yo = float(yo)
        a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
        b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
        c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
        g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                               + c*((y-yo)**2)))
        return g.ravel()

# Get image:
d = pyfits.getdata('nestor/results/master_ao_no_flats_400_hs.fits')
# Set scale:
scale = 0.016 # arcsec/pix

# Measure centroid and FWHM of the central star:
A,x0,y0,sigma_x,sigma_y,theta,offset = Utils.get_centroid(d,mode = 'non-subimg')

# Define magnitude difference:
DeltaMag = 5
# Define distance (in pixels) and angle at which you want to insert the fake source:
r = 130.
theta = 160.0

Astar = A*10**(-DeltaMag/2.51)

print 'sigma_x: ',sigma_x*scale,'arcsec','sigma_y:',sigma_y*scale, 'arcsec'
print 'fwhm_x: ',sigma_x*scale*2.355,'arcsec','fwhm_y:',sigma_y*scale*2.355, 'arcsec'
# Plot best-fit:
x = np.arange(d.shape[0])
y = np.arange(d.shape[1])
x,y = np.meshgrid(x,y)

data = twoD_Gaussian((x,y),Astar,x0+r*np.cos(theta*np.pi/180.),y0+r*np.sin(theta*np.pi/180.),sigma_x,sigma_y,theta,0.0)
ndata = data.reshape(d.shape[0],d.shape[1])
pyfits.PrimaryHDU(ndata+d).writeto('fit.fits')


