# download data from pre-loaded github folder
from shutil import rmtree
from os import path
    
dir_images = 'WEW2024_database/LASCO_data/'

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sunpy.map
import glob  # To fetch all files from a directory

map_list = sorted(glob.glob(dir_images+'*.fts'))

lasco_map = sunpy.map.Map(map_list[0])

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=lasco_map)
lasco_map.plot(axes=ax, cmap = plt.get_cmap('soholasco3'))
lasco_map.draw_limb()
plt.show()

#get the time of the map
print(lasco_map.date)

from os import makedirs as newdir
dir_png = 'My_pics_LASCO'
newdir(dir_png, exist_ok = True)

for index, picture in enumerate(map_list):
  lasco_map = sunpy.map.Map(picture)
  lasco_map.plot()
  plt.savefig(f'{dir_png}/{index:03d}.png', dpi=150)
  plt.clf()
  
  print(f'Saved figure {index}')

  import cv2
img_array = []

for filename in sorted(glob.glob(f'{dir_png}/*.png')):
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width,height)
    img_array.append(img)

out = cv2.VideoWriter('solar_vid_coronograph.mp4',cv2.VideoWriter_fourcc(*'mp4v'), 5, size)

for i in range(len(img_array)):
    out.write(img_array[i])

out.release()