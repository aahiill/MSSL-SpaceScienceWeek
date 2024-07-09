from astropy import units as u
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from ccdproc import CCDData, combine, ccd_process, subtract_bias, subtract_dark
from ccdproc.combiner import Combiner
from astroscrappy import detect_cosmics
import numpy as np

def load_ccd_data(file_path):
    return CCDData.read(file_path, unit='adu')

obj = input("Enter which object (1, 2, 3)   :   ")

f1 = 'C:/Users/aahil/Desktop/Code/MSSL_SpaceScienceWeek/Astro_Physics/Object'+obj+'/file1.fits'
f2 = 'C:/Users/aahil/Desktop/Code/MSSL_SpaceScienceWeek/Astro_Physics/Object'+obj+'/file2.fits'
f3 = 'C:/Users/aahil/Desktop/Code/MSSL_SpaceScienceWeek/Astro_Physics/Object'+obj+'/file3.fits'
f4 = 'C:/Users/aahil/Desktop/Code/MSSL_SpaceScienceWeek/Astro_Physics/Object'+obj+'/file4.fits'
f5 = 'C:/Users/aahil/Desktop/Code/MSSL_SpaceScienceWeek/Astro_Physics/Object'+obj+'/file5.fits'
f6 = 'C:/Users/aahil/Desktop/Code/MSSL_SpaceScienceWeek/Astro_Physics/Object'+obj+'/file6.fits'
f7 = 'C:/Users/aahil/Desktop/Code/MSSL_SpaceScienceWeek/Astro_Physics/Object'+obj+'/file7.fits'
f8 = 'C:/Users/aahil/Desktop/Code/MSSL_SpaceScienceWeek/Astro_Physics/Object'+obj+'/file8.fits'

f1_data = load_ccd_data(f1)
f2_data = load_ccd_data(f2)
f3_data = load_ccd_data(f3)
f4_data = load_ccd_data(f4)
f5_data = load_ccd_data(f5)
f6_data = load_ccd_data(f6)
f7_data = load_ccd_data(f7)
f8_data = load_ccd_data(f8)


plt.imshow(np.log10(f1_data), cmap = 'gray')
plt.title('f1')
plt.show()
plt.imshow(np.log10(f2_data), cmap = 'gray')
plt.title('f2')
plt.show()
plt.imshow(np.log10(f3_data), cmap = 'gray')
plt.title('f3')
plt.show()
plt.imshow(np.log10(f4_data), cmap = 'gray')
plt.title('f4')
plt.show()
plt.imshow(np.log10(f5_data), cmap = 'gray')
plt.title('f5')
plt.show()
plt.imshow(np.log10(f6_data), cmap = 'gray')
plt.title('f6')
plt.show()
plt.imshow(np.log10(f7_data), cmap = 'gray')
plt.title('f7')
plt.show()
plt.imshow(np.log10(f8_data), cmap = 'gray')
plt.title('f8')
plt.show()

blue = f1 # this is the science file for the blue (B#639) filter.
visible = f2 # this is the science file for the visible (V#641) filter.
red = f3 # this is the science file for the red (R#642) filter.
blue_flat = f4 # this is the flat file for the blue (B#639) filter.
visible_flat = f5 # this is the flat file for the visible (V#641) filter.
red_flat = f6 # this is the flat file for the red (R#642) filter.
bias = f7 # this is the bias file.
dark = f8 # this is the dark file.

def process_filter(flat, science, bias, dark):
    """
    Reads the science, flat, biad and dark files for a specific filter and
    performs the calibration procedure

    Parameters
    ----------
    flat : string
        Path to flat file.
    science : string
        path to science file.
    bias : string
        path to bias file.
    dark : string
        path to dark file.


    Returns
    -------
    flat_corrected : CCDData class
        the calibrated data

    """
    #Load the data
    flat_frames = load_ccd_data(flat)
    science_data = load_ccd_data(science)
    bias_frames = load_ccd_data(bias)
    dark_frames = load_ccd_data(dark)

    #Obtain the master dark calibration data
    bias_subtracted_dark = subtract_bias(dark_frames, bias_frames)
    master_dark = bias_subtracted_dark

    #Correct for bias
    bias_corrected = subtract_bias(science_data, bias_frames)
    #Correct for dark
    dark_corrected = subtract_dark(bias_corrected, master_dark, dark_exposure=master_dark.meta['EXPTIME'] * u.second, data_exposure=science_data.meta['EXPTIME'] * u.second)
    #Correct for flat
    flat_corrected = ccd_process(dark_corrected, master_flat = flat_frames)
    flat_corrected.write('processes_file.fits', overwrite = True)

    data = fits.getdata('processes_file.fits')
    #Remove comsic rays and bad pixels
    cr_mask, clean_data = detect_cosmics(data, sigclip=4.5, sigfrac = 0.3, objlim=5.0)

    return flat_corrected

#Calibrate the data for all filters
corrected_blue = process_filter(blue_flat, blue, bias, dark)
corrected_red = process_filter(red_flat, red, bias, dark)
corrected_visible = process_filter(visible_flat, visible, bias, dark)

#Transform any infinities into np.NaNs
corrected_blue.data[np.isinf(corrected_blue.data)] = np.nan
corrected_red.data[np.isinf(corrected_red.data)] = np.nan
corrected_visible.data[np.isinf(corrected_visible.data)] = np.nan

#Plot the corrected data
plt.imshow(np.log10(corrected_blue), cmap = 'Blues')
plt.title('Blue')
plt.show()
plt.imshow(np.log10(corrected_red), cmap = 'Greens')
plt.title('Visible')
plt.show()
plt.imshow(np.log10(corrected_visible), cmap = 'Reds')
plt.title('Red')
plt.show()

#Plot a histogram for each filter to assess the photon distribution across the image.
plt.hist(np.ravel(corrected_blue.data), bins = 100, range = (0, 40e2), color = 'blue')
plt.title('Blue')
plt.show()
plt.hist(np.ravel(corrected_visible.data), bins = 100, range = (0, 40e2), color = 'green')
plt.title('Visible')
plt.show()
plt.hist(np.ravel(corrected_red.data), bins = 100, range = (0, 40e2), color = 'red')
plt.title('Red')
plt.show()

def scale_image(image, min_z, max_z, log_scale=False):

    image = np.clip(image, min_z, max_z)

    if log_scale:
      image_log = np.log10(image)
      min_z_log = np.log10(min_z)
      max_z_log = np.log10(max_z)

      image_scaled = (image_log - min_z_log) / (max_z_log - min_z_log)

    else:
      image_scaled = (image - min_z) / (max_z - min_z)

    return image_scaled

#Change these parameters
blue_min = 1
blue_max = 1000
visible_min = 1
visible_max = 1000
red_min = 1
red_max = 1000
blue_ratio = 1
visible_ratio = 1
red_ratio = 1
use_log = True

#Scale each filter
blue_scaled = blue_ratio * scale_image(corrected_blue, blue_min, blue_max, log_scale = use_log)
visible_scaled = visible_ratio * scale_image(corrected_visible, visible_min, visible_max, log_scale = use_log)
red_scaled = red_ratio * scale_image(corrected_red, red_min, red_max, log_scale = use_log)

#Plot the scaled images for each filter
plt.imshow(blue_scaled, cmap = 'Blues')
plt.show()
plt.imshow(visible_scaled, cmap = 'Greens')
plt.show()
plt.imshow(red_scaled, cmap = 'Reds')
plt.show()

#Combine the filters into an rgb_image
rgb_image = np.stack((red_scaled, visible_scaled, blue_scaled), axis = -1)
rgb_image = np.clip(rgb_image, 0, 1)

#Plot the rgb image
plt.imshow(rgb_image)
plt.show()