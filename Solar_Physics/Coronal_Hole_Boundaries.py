# Coronal Hole Boundaries

import glob

from os import makedirs as newdir
from os import path

import sunpy.map
from sunpy.net import attrs as a
from sunpy.net import Fido

from sunpy.map.maputils import all_pixel_indices_from_map,sample_at_coords,all_coordinates_from_map, coordinate_is_on_solar_disk

from aiapy.calibrate import register, update_pointing, correct_degradation

import numpy as np

from scipy import ndimage
from scipy.stats import mode

import matplotlib.pyplot as plt

import astropy.time
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord

import warnings
warnings.simplefilter('ignore')

# download data from pre-loaded github folder
from shutil import rmtree
from os import path

dir_images = 'WEW2024_database/Sampledata_AIA193/'

filelist_AIA = sorted(glob.glob(dir_images+'*.fits'))
# Select first file (first element => index 0) 
# You can change the image file by changing the index number
fileAIA = filelist_AIA[0]
aiamap = sunpy.map.Map(fileAIA)

#Function use to prepare map obtained from AIA
def aiaprep(aiamap):
    aiamap_updated_pointing = update_pointing(aiamap) #update metadata to recent pointing
    aiamap_registered = register(aiamap_updated_pointing) #registered: rotate align w/ solar north and rescale to match 0.6 pixel resolution
    aiamap_norm = aiamap_registered/aiamap_registered.exposure_time #Normalized exposure time for each image: unit (DN/pixel/s)
    aiamap_corr = correct_degradation(aiamap_norm) #Correct Degradation effect of Instrument
    aiamap_rot = aiamap_corr.rotate()

    return aiamap_rot

aiamap_prep = aiaprep(aiamap) #input is full sun image file

aiamap_prep.peek()

# Choose submap region
bottom_left_x = -1000 #x coordinate of bottomleft corner
bottom_left_y = -800 #y coordinate of bottomleft corner
bottom_left = SkyCoord(bottom_left_x * u.arcsec, bottom_left_y* u.arcsec, frame= aiamap_prep.coordinate_frame)
# Define width and height for our zoomed-in map
w = 800 #Width
h = 1000 #Height
# Crop our map using submap function
subaiamap = aiamap_prep.submap(bottom_left, width = w*u.arcsec, height = h*u.arcsec)
subaiamap_mask = aiamap_prep.submap(bottom_left, width = w*u.arcsec, height = h*u.arcsec)

# Define the number of bins: The number of bins will affect the shape of distribution
num_bins = 1000
bins = np.linspace(subaiamap.min(), subaiamap.max(), num_bins)

# Plot Histogram
fig = plt.figure()
plt.hist(subaiamap.data.ravel(), bins=bins, label='Histogram', histtype='step' ,density = True)
plt.xlabel('Intensity')
plt.axvline(subaiamap.data.mean(),
            label='mean={:.2f}'.format(subaiamap.data.mean()), color='green')
plt.axvline(np.median(subaiamap.data),
            label='median={:.2f}'.format(np.median(subaiamap.data)), color='Orange')
plt.legend(loc=9)
# Set the limit of x axis to zoomed in at specific range 
# (Hint: appropriate range is x-axis can be helpful)
plt.xlim(0,200)
plt.show()

# Change the value of thr to your guess
thr = input("Enter the intensity threshold value (40): ") 
print(thr)

#Function used to mask the area outside solar disc
def mask_solardisc(subaiamap):
    hpc = all_coordinates_from_map(subaiamap)
    mask = coordinate_is_on_solar_disk(hpc)
    return mask

#Define Coronal Hole Boundary - Return results as the array where CH pixels are defined as 1 and anywhere else are defined as 0
def define_CHB (aiamap, subaiamap, subaiamap_mask, thr, filt_order=5): 
    mask_disk = mask_solardisc(subaiamap)
    mask_disk_full = mask_solardisc(aiamap)
    #Make the data outside solar disk = 0
    aiamap_data_disk = np.where(mask_disk_full, aiamap.data, 0)
    #Select CHB according to given threshold
    mask = subaiamap.data > thr
    subaiamap_mask.mask = mask
    #Gaussian Filtering to Connect Smaller Region and smooth boundary
    data1 = ndimage.gaussian_filter(subaiamap.data * ~mask, filt_order) #gaussian filtering
    #labels the area define in CH
    subaiamap_new = sunpy.map.Map(data1, subaiamap.meta)
    labels, n = ndimage.label(subaiamap_new.data)
    #Select only biggest CH
    vallabel = mode(labels[np.nonzero(labels)])
    labels_true = np.where(labels == vallabel[0], labels, 0)
    #Fill Bright spots in CH (We neglected them)
    labels_true = ndimage.binary_fill_holes(labels_true).astype(int)
    # Check only coordinate on disk
    labels_true = np.where(mask_disk, labels_true,0)
    return labels_true

# input are (full sun map, submap, a copy of submap(mask), threshold value)
CHB = define_CHB(aiamap, subaiamap, subaiamap_mask, thr=40)

fig = plt.figure()
ax = fig.add_subplot(1,1,1,projection = subaiamap)
subaiamap.plot()
# ax.contour will draw a contour around boundary between 0 and 1 in CHB array, which is the boundary of CH
ax.contour(CHB, colors='black', linewidths = 0.5) #You can change color and linewidths

dir_png = 'My_pics_CHB'
newdir(dir_png, exist_ok = True)
# Save plots using savefig we can also have option such as choosing name, select the resolution or eliminate white space
fig.savefig(dir_png+'/CHB_AIA193.png', dpi=200, bbox_inches='tight')

def getboundarycoord(labels, subaia):
    maskinner =  ndimage.binary_erosion(labels.tolist())
    labels[maskinner] = 0
    CHB_pixel = np.nonzero(labels)
    CHB_pixel_array = np.column_stack((CHB_pixel[1],CHB_pixel[0]))
    xpixelaia = CHB_pixel_array[:,0]
    ypixelaia = CHB_pixel_array[:,1]
    subaiamapwcs = subaia.wcs
    aia_subarea_sky = pixel_to_skycoord(xpixelaia,ypixelaia,subaiamapwcs,origin=0)
    return aia_subarea_sky

CHB_coord = getboundarycoord(CHB, subaiamap)

fig = plt.figure()
ax = fig.add_subplot(1,1,1,projection = aiamap)
aiamap.plot()
# ax.plot_coord will plot the coordinate, in this case is our CH boundary
ax.plot_coord(CHB_coord, '.', color='black', markersize = 0.1) #You can change color and linewidths
plt.show()

