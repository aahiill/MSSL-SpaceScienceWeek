from sunpy.net import Fido, attrs as a
from astropy import units as u  

from datetime import datetime   # Python's small Time module

import warnings
warnings.simplefilter('ignore')  # We will ignore warnings for now

valid_wavelengths = [94, 131, 171, 193, 211, 304, 335, 1600, 1700]  # Note that 16/1700 are white light and not as useful here!


Inst = 'aia'  # Here we will use the Atmospheric Imaging Assembly: aia
# It is onboard the Solar Dynamics Observatory. 
#Provides us with images of the Sun almost 24-7!


# Now choose the period we would like to get our data from. Please don't do too much at once!

start_time = datetime(2012,8,31,18,0)  # Year,Month,Day,Hour,Minute
end_time = datetime(2012,8,31,23,59)

time_range = a.Time(start_time, end=end_time)

lambda_ang = valid_wavelengths[2] * u.Angstrom  # Choose the 171 Angstrom wavelength as an example (Extreme Ultraviolet)


samplerate_hours = 12 * u.minute  # We need to give our images a cadence for collection. 
# Please keep in mind that aia takes an image in all of its wavelengths every 12 seconds.
# Choose a cadence depending on what you want to bring out!

aia_query = time_range & a.Instrument.aia & a.Wavelength(lambda_ang) & a.Sample(samplerate_hours)

fido = Fido.search(aia_query)

print(fido)

# download data from pre-loaded github folder [git clone https://github.com/nawinnova/WEW2024_database.git]
from shutil import rmtree
from os import path

wavelength = str(input("Enter database wavelength (171, 193)  :   "))

    
dir_images = 'WEW2024_database/Sampledata_AIA'+wavelength+'/'

import matplotlib.pyplot as plt
import sunpy.map
import glob  # To fetch all files from a directory
from os import makedirs as newdir


map_list = sorted(glob.glob(dir_images + '*.fits'))
 
for picture in map_list[0:3]:
  
  aia_map = sunpy.map.Map(picture)
  aia_map.peek()
  
input('Happy with the results? Press enter again in this box to save all of the figures as png images')


dir_png = 'My_pics_AIA'+wavelength
newdir(dir_png, exist_ok = True)

for index, picture in enumerate(map_list):
  aia_map = sunpy.map.Map(picture)
  aia_map.plot()
  plt.savefig(f'{dir_png}/{index:03d}.png', dpi=150)
  plt.clf()
  
  print(f'Saved figure {index}')

  import cv2  # Here, we import OpenCV, a powerful image and video processing tool
from os import getcwd


img_array = []

for filename in sorted(glob.glob(f'{dir_png}/*.png')):
  img = cv2.imread(filename)
  height, width, layers = img.shape
  size = (width,height)
  img_array.append(img)

out = cv2.VideoWriter('solar_vid_disk_AIA'+wavelength+'.mp4',cv2.VideoWriter_fourcc(*'MP4V'), 6, size)
 
for i in range(len(img_array)):
    out.write(img_array[i])
out.release()
