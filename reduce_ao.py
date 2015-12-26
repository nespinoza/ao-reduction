import os
import numpy as np
import matplotlib.pyplot as plt
import Utils
import glob
import pyfits
from scipy.ndimage.interpolation import shift as img_shift
from scipy.ndimage.filters import median_filter
plt.style.use('ggplot')

############################ Parameter #######################
reduction_name = 'nestor'
scale = 0.016 # arcsec/pix
foldername = '/Volumes/SeagateEHD/data/KEPLER/AO/Nestor_star/'
#foldername = '/Volumes/SeagateEHD/data/AO/clio_20151206_07/'
half_size = 400

use_sky_flats = False
median_filtering = True
window = 200
# Get the images:
#image_fnames = glob.glob(foldername+'target2_south_00*.fit')#'Nestor1_00*.fit')
image_fnames = glob.glob(foldername+'Nestor1_00*.fit')
# Get the darks:
#dark_frames = glob.glob(foldername+'dark00*.fit')
dark_frames = glob.glob(foldername+'darkN20*.fit')
# Get sky flats:
sky_frames = glob.glob(foldername+'skyflat00*.fit')

###############################################################
print '\t > Pre-processing...'
if not os.path.exists(reduction_name):
   os.mkdir(reduction_name)

if use_sky_flats:
    calib_folder = reduction_name+'/calibration_w_flats'
else:
    calib_folder = reduction_name+'/calibration_wo_flats'

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
    bad_pixel_map = pyfits.getdata(calib_folder+'/badpix_fullframe.fit')
    
# Get bad pixels to zero, good ones to one:
bad_pixel_map = np.abs(bad_pixel_map-1)

out_images_folder = calib_folder+'/corrected_images/'
if not os.path.exists(out_images_folder):
    os.mkdir(out_images_folder)

print '\t > Calibrating each image...'
for i in range(len(image_fnames)):
    data,h = pyfits.getdata(image_fnames[i],header=True)
    data = ((data-median_dark)/median_sky)*bad_pixel_map
    data = data - np.median(data)
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

# First median substract all the images and get their centroids:
print '\t > Calculating centroids...'
x_centroids = np.zeros(len(all_images))
y_centroids = np.zeros(len(all_images))
for i in range(len(all_images)):
    # Substract the sky:
    all_images[i] = all_images[i]-median_image
    # If median filtering is activated, generate it:
    if median_filtering:
        mf = median_filter(all_images[i],size=window)
        all_images[i] = all_images[i] - mf
    # Get centroid of brightest star on sky-substracted image:
    x_centroids[i],y_centroids[i] = Utils.get_centroid(all_images[i])

# Cluster images with similar centroids together (they have the same nodding position):
print '\t > Clustering nodding patterns...'
patterns = []
all_idx = []
init_i = 0
while True:
    c_idx = []
    c_x,c_y = x_centroids[init_i],y_centroids[init_i]
    for i in range(len(all_images)):
        if i not in all_idx:
            dist = np.sqrt((c_x-x_centroids[i])**2 + (c_y-y_centroids[i])**2)
            if dist<3.:
                c_idx.append(i)
                all_idx.append(i)
            else:
                init_i = i
    patterns.append(c_idx)
    if len(all_idx) == len(all_images):
       break

if not os.path.exists(reduction_name+'/results'):
    os.mkdir(reduction_name+'/results')
# Median combine each nodding position separately. First, get rotated images in 
# a common reference frame:
print '\t Rotating, shifting and combining...'
rotated_frames = len(all_images)*[[]]
for i in range(len(all_images)):
    diff_image = all_images[i]
    x,y = x_centroids[i],y_centroids[i]
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
    pyfits.PrimaryHDU(shifted_image[center_x-half_size:center_x+half_size,\
                      center_y-half_size:center_y+half_size]).writeto(\
                      reduction_name+'/results/image'+str(i)+'.fits')

    rotated_frames[i] = shifted_image[center_x-half_size:center_x+half_size,\
                                      center_y-half_size:center_y+half_size]

for i in range(len(patterns)):
    pattern = patterns[i]
    for j in range(len(pattern)):
        if j == 0:
           all_c_data = np.copy(rotated_frames[pattern[j]])
        else:
           all_c_data = np.dstack((all_c_data,np.copy(rotated_frames[pattern[j]])))

    pyfits.PrimaryHDU(np.median(all_c_data,axis=2)).writeto(reduction_name+'/results/pattern'+str(i)+'.fits')

    if i == 0:
       master_images = np.median(all_c_data,axis=2)
    else:
       master_images = np.dstack((master_images,np.median(all_c_data,axis=2)))

#for i in range(len(patterns)):
#    if i == 0:
#        master_image = np.copy(shifted_image[center_x-half_size:center_x+half_size,\
#                                             center_y-half_size:center_y+half_size])
#        master_count = np.ones(master_image.shape)
#    else:
#        master_image = master_image + shifted_image[center_x-half_size:center_x+half_size,\
#                                                    center_y-half_size:center_y+half_size]
#        add_matrix = np.zeros(master_image.shape)
#        xx,yy = np.where(shifted_image[center_x-half_size:center_x+half_size,\
#                                      center_y-half_size:center_y+half_size]>0.0)
#        add_matrix[xx,yy] = np.ones(len(xx))
#        master_count = master_count + add_matrix
#
#im = plt.imshow(master_image/master_count)
result = np.median(master_images,axis=2)
im = plt.imshow(result)
im.set_clim(-50,50)
plt.show()

if not os.path.exists(reduction_name+'/results'):
   os.mkdir(reduction_name+'/results')

if use_sky_flats:
    pyfits.PrimaryHDU(result).writeto(reduction_name+'/results/master_ao_w_flats.fits')
else:
    pyfits.PrimaryHDU(result).writeto(reduction_name+'/results/master_ao_no_flats.fits')
