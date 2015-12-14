import scipy.optimize as opt
import numpy as np

def get_centroid(image, A_init, x_0_init, y_0_init, sigma_x_init, sigma_y_init, theta_init, offset_init):
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

    # Get x and y of image:
    x = np.arange(image.shape[0])
    y = np.arange(image.shape[1])
    x, y = np.meshgrid(x, y)

    # Define initial guess:
    initial_guess = (A_init, x_0_init, y_0_init, sigma_x_init, sigma_y_init, theta_init, offset_init)

    flattened_image = image.flatten()
    popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), flattened_image, p0=initial_guess)

    return popt[1],popt[2]
   

def get_init_params(image): 


