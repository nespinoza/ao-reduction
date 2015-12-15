import os
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
half_size = 250

use_sky_flats = False
# Get the images:
image_fnames = glob.glob(foldername+'Nestor1_00*.fit')
# Get the darks:
dark_frames = glob.glob(foldername+'darkN200*.fit')
# Get sky flats:
sky_frames = glob.glob(foldername+'skyflat00*.fit')

###############################################################

if use_sky_flats:
    calib_folder = 'calibration_w_flats'
else:
    calib_folder = 'calibration_wo_flats'

if not os.path.exists(calib_folder):
    os.mkdir(calib_folder)

# Create master dark:
for i in range(len(dark_frames)):
    data,h = pyfits.getdata(dark_frames[i],header=True)
    if i == 0:
        dark_images = np.copy(data)
    else:
        dark_images = np.dstack((dark_images,np.copy(data)))

median_dark = np.median(dark_images,axis=2)

if use_sky_flats:
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
else:
    median_sky = 1.0

# Now create median image from all the science frames:
all_images = [] 
all_headers = []

# Get bad pixel map:
try:
    bad_pixel_map = pyfits.getdata(calib_folder+'/badpix_fullframe.fit')
except:
    print 'Bad pixel map not found, downloading it...'
    os.system('wget http://zero.as.arizona.edu/groups/clio2usermanual/wiki/6d927/attachments/1242c/badpix_fullframe.fit')
    os.system('mv badpix_fullframe.fit '+calib_folder+'/badpix_fullframe.fit')
    bad_pixel_map = pyfits.getdata(calib_folder+'badpix_fullframe.fit')
    
# Get bad pixels to zero, good ones to one:
bad_pixel_map = np.abs(bad_pixel_map-1)

out_images_folder = calib_folder+'/corrected_images/'
if not os.path.exists(out_images_folder):
    os.mkdir(out_images_folder)

for i in range(len(image_fnames)):
    data,h = pyfits.getdata(image_fnames[i],header=True)
    data = ((data-median_dark)/median_sky)*bad_pixel_map
    if not os.path.exists(out_images_folder+image_fnames[i].split('/')[-1]):
        pyfits.PrimaryHDU(data).writeto(out_images_folder+image_fnames[i].split('/')[-1])
    if i == 0:
        all_data = np.copy(data)
    else:
        all_data = np.dstack((all_data,np.copy(data)))
    all_images.append(np.copy(data))
    all_headers.append(h)

try:
    median_image = pyfits.getdata(out_images_folder+'median_image.fits')
except:
    median_image = np.median(all_data,axis=2)
    pyfits.PrimaryHDU(median_image).writeto(out_images_folder+'median_image.fits')

# Now get final image by rotating and shifting the images. The 
# rotation is done according to what is posted here: 
# http://zero.as.arizona.edu/groups/clio2usermanual/wiki/6d927/Calibration_Data.html
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
    center_x = int(shifted_image.shape[0]/2.)
    center_y = int(shifted_image.shape[1]/2.)
    if i == 0:
        master_image = np.copy(shifted_image[center_x-half_size:center_x+half_size,\
                                             center_y-half_size:center_y+half_size])
        master_count = np.ones(master_image.shape)
    else:
        master_image = master_image + shifted_image[center_x-half_size:center_x+half_size,\
                                                    center_y-half_size:center_y+half_size]
        add_matrix = np.zeros(master_image.shape)
        xx,yy = np.where(shifted_image[center_x-half_size:center_x+half_size,\
                                      center_y-half_size:center_y+half_size]>0.0)
        add_matrix[xx,yy] = np.ones(len(xx))
        master_count = master_count + add_matrix
    #im = plt.imshow(shifted_image)
    #im.set_clim(-50,50)
    #plt.plot(xoff,yoff,'ro')
    #plt.show()

im = plt.imshow(master_image/master_count)
im.set_clim(-50,50)
plt.show()
if use_sky_flats:
    pyfits.PrimaryHDU(master_image/master_count).writeto('master_ao_w_flats.fits')
else:
    pyfits.PrimaryHDU(master_image/master_count).writeto('master_ao_no_flats.fits')
