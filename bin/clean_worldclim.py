import rasterio
import numpy as np

# Load a specific WorldClim GeoTIFF
dataset = rasterio.open('path/to/wc2.1_30s_prec_01.tif')
precipitation_january = dataset.read(1)

