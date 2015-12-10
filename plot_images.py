import numpy as np
import matplotlib.pyplot as plt
import glob
import pyfits
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
sky_frames = glob.glob(foldername+'skyflat200*.fit')

# Create master dark:
for i in range(len(dark_frames)):
    data,h = pyfits.getdata(dark_frames[i],header=True)
    if i == 0:
        dark_images = np.copy(data)
    else:
        dark_images.dstack((dark_images,np.copy(data)))

median_dark = np.median(dark_image,axis=2)

# Create master sky:
for i in range(len(dark_frames)):
    data,h = pyfits.getdata(sky_frames[i],header=True)
    data = data - median_dark
    if i == 0:
        sky_images = np.copy(data/np.median(data))
    else:
        sky_images.dstack((sky_images,np.copy(data/np.median(data))))

median_sky = np.median(sky_image,axis=2)

all_images = [] 
for i in range(len(image_fnames)):
    data,h = pyfits.getdata(image_fnames[i],header=True)
    if i == 0:
        all_data = (data-median_dark)/median_sky
    else:
        all_data = np.dstack((all_data,(data-median_dark)/median_sky))
    all_images.append(np.copy(data))

median_image = np.median(all_data,axis=2)
for i in range(len(all_images)):
    diff_image = all_images[i]-median_image
    im = plt.imshow(diff_image)
    sigma = np.sqrt(np.var(diff_image))
    im.set_clim(np.median(diff_image)-sigma,np.median(diff_image)+sigma)
    x,y = np.where(diff_image == np.max(diff_image))
    plt.plot([y,y],[x,x-1./scale],'-')
    plt.plot([y,y],[x,x-4./scale],'--')
    plt.plot([y,y],[x,x+1./scale],'-')
    plt.plot([y,y],[x,x+4./scale],'--')
    plt.plot([y,y-1./scale],[x,x],'-')
    plt.plot([y,y-4./scale],[x,x],'--')
    plt.plot([y,y+1./scale],[x,x],'-')
    plt.plot([y,y+4./scale],[x,x],'--')
    plt.show()

