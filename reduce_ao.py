import numpy as np
import matplotlib.pyplot as plt
import Utils
import glob
import pyfits
from scipy.ndimage.interpolation import shift as img_shift
plt.style.use('ggplot')

############################ Parameter #######################
scale = 0.016 # arcsec/pix
foldername = '/Volumes/SeagateEHD/data/KEPLER/AO/Nestor_star/'

##############################################################
# Get the images:
image_fnames = glob.glob(foldername+'Nestor1_00*.fit')
# Get the darks:
dark_frames = glob.glob(foldername+'darkN200*.fit')
# Get sky flats:
sky_frames = glob.glob(foldername+'skyflat00*.fit')

# Create master dark:
for i in range(len(dark_frames)):
    data,h = pyfits.getdata(dark_frames[i],header=True)
    if i == 0:
        dark_images = np.copy(data)
    else:
        dark_images = np.dstack((dark_images,np.copy(data)))

median_dark = np.median(dark_images,axis=2)

# Create master sky:
for i in range(len(sky_frames)):
    data,h = pyfits.getdata(sky_frames[i],header=True)
    data = data - median_dark
    if i == 0:
        sky_images = np.copy(data/np.median(data))
    else:
        sky_images = np.dstack((sky_images,np.copy(data/np.median(data))))

# Correct for zeros on sky image:
median_sky = np.median(sky_images,axis=2)
x_z,y_z = np.where(median_sky == 0)
for i in range(len(x_z)):
    median_sky[x_z,y_z] = 1.0

# Now create median image from all the science frames:
all_images = [] 
all_headers = []
for i in range(len(image_fnames)):
    data,h = pyfits.getdata(image_fnames[i],header=True)
    data = (data-median_dark)/median_sky
    if i == 0:
        all_data = np.copy(data)
    else:
        all_data = np.dstack((all_data,np.copy(data)))
    all_images.append(np.copy(data))
    all_headers.append(h)

median_image = np.median(all_data,axis=2)

# Graph every image:
for i in range(len(all_images)):
    # Substract the sky:
    diff_image = all_images[i]-median_image
    # Get centroid of brightest star on sky-substracted image:
    x,y = Utils.get_centroid(diff_image)
    # Rotate the image:
    rotated_image = Utils.rotate_image(diff_image,all_headers[i])
    # Get the centers of the original image:
    y_c = diff_image.shape[0]/2.
    x_c = diff_image.shape[1]/2.
    # Get centers of the rotated image (everything is rotated w/r to this):
    yoff = rotated_image.shape[0]/2.
    xoff = rotated_image.shape[1]/2.
    # Calculate the centroid of the rotated star:
    xrot,yrot = Utils.rotate_point(y+(xoff-x_c),x+(yoff-y_c),all_headers[i],xoff,yoff)
    # Shift the image in order to leave the centroid at the center of the image:
    shifted_image = img_shift(rotated_image,[(yoff-yrot),(xoff-xrot)])
    im = plt.imshow(shifted_image)
    im.set_clim(-50,50)
    plt.plot(xoff,yoff,'ro')
    plt.show()

