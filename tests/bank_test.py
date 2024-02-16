#%%
import os
from rastertools import delineate_bank_lines
from rastertools import smooth_line_geom

dem_file = '../data/brazos/input/dem_clip.tif'
banks = delineate_bank_lines(dem_file)
banks.plot()
