# -*- coding: utf-8 -*-
import scipy.optimize as opt
import numpy as np

def get_centroid(original_image):
    """
    get_centroid
    ------------

    This algorithm assumes your image has only one bright source. 
    It calculates the centroid of such source. Idea taken from here: 
    http://stackoverflow.com/questions/21566379/fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m
    """
    def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
        xo = float(xo)
        yo = float(yo)    
        a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
        b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
        c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
        g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                               + c*((y-yo)**2)))
        return g.ravel()

    # Get initial guess for centroids and subimg:
    x_0_init,y_0_init,x_init_subimg,y_init_subimg,image = get_init_centroid(original_image)
    sigma_x_init = 10.
    sigma_y_init = 10.
    theta_init = 0.0
    offset_init = 0.0

    # Get x and y of image:
    x = np.arange(image.shape[0])
    y = np.arange(image.shape[1])
    x, y = np.meshgrid(x, y)

    A_init = np.sum(image)/(2.*np.pi*sigma_x_init*sigma_y_init)

    # Define initial guess:
    initial_guess = (A_init, x_0_init, y_0_init, sigma_x_init, sigma_y_init, theta_init, offset_init)

    flattened_image = image.flatten()
    popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), flattened_image, p0=initial_guess)

    return x_init_subimg+popt[1],y_init_subimg+popt[2]


from scipy.ndimage.filters import median_filter
def get_init_centroid(image,half_size = 100,median_window = 30): 
    """
    Initial guess for the centroid. The idea is pretty simple: 
    get rid of deviating pixels with a median filter; hence, create 
    a smooth version of the image with no outliers. Then, get as first 
    guess of the centroid the maximum count value (i.e., the brightest 
    star in the image).

    The idea is then to cut a sub-image around this star, and fit a gaussian 
    to that.
    """
    mf = median_filter(image,size=median_window)
    x0,y0 = np.where(mf == np.max(mf))
    x0,y0 = x0[0],y0[0]
    x_init = np.max([0,int(x0)-half_size])
    x_end = np.min([mf.shape[0],int(x0)+half_size])
    y_init = np.max([0,int(y0)-half_size])
    y_end = np.min([mf.shape[1],int(y0)+half_size])

    # Get subimg:
    sub_img = image[x_init:x_end,\
                    y_init:y_end]

    return x0-x_init,y0-y_init,x_init,y_init,sub_img

# Rotate the the image according to DEROT = ROTOFF - 180 + NORTHCLIO 
# (http://zero.as.arizona.edu/groups/clio2usermanual/wiki/6d927/Calibration_Data.html)
from scipy.ndimage.interpolation import rotate
NorthClio = -1.80
def rotate_image(image,header):
    angle = header['ROTOFF'] - 180. + NorthClio
    return rotate(image,header['ROTOFF'] - 180. + NorthClio,reshape=True)

# Given a rotation about x_off and y_off, rotate a point in the same angle 
# as the image given its x,y coordinates.
def rotate_point(x,y,header,x_off,y_off):
    angle = -(header['ROTOFF'] - 180. + NorthClio)*np.pi/180.
    return ((x-x_off)*np.cos(angle) - (y-y_off)*np.sin(angle)) + x_off,\
           ((x-x_off)*np.sin(angle) + (y-y_off)*np.cos(angle)) + y_off
