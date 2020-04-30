
# This scriprt loop all GeoTIFF files in a folder with a specific fligh line name, assuming they are all trait files.
# Then convert the traits files to BSQ-interleaved GeoTIFF. Meanwhile, it tries to detect the Nodata value by probing all four corners of the image. 

# Usage: python remake_trait_geotiff.py  input_dir output_dir   flight_line_name
# e.g. "python remake_trait_geotiff.py  /hyspiri/2016/in/ /hyspiri/2016/out/   f160616t01p00r17"

import sys,  os, glob
import gdal, osr
import numpy as np

from collections import Counter

from gdalconst import *


out_nodata_val = -9999
out_suffix = "_v2.tif" # output suffix of the geotiff

inpath = sys.argv[1]
outpath = sys.argv[2]
flight_name = sys.argv[3]

pat = os.path.join(inpath,flight_name+"*.tif")


# Read file and save it in memory
# It looks for NoData based on the values at the four corners
def read_img(srcfile):
  
  src_ds=gdal.Open(srcfile)
  src_wkt=src_ds.GetProjection() #wkt
  src_geotransform =src_ds.GetGeoTransform()

  #src_shape = (src_ds.RasterYSize, src_ds.RasterXSize)
  
  # assuming there are two bands, trait model mean and model stdev
  trait_avg = src_ds.GetRasterBand(1).ReadAsArray()
  trait_std = src_ds.GetRasterBand(2).ReadAsArray()
  
  # Initialized with large enough values
  trait_avg_nodata=9999
  trait_std_nodata=9999
  

  # Get the values for the four corners, and try to find the mode of these four, then set the mode value as the final NoData value
  c1=Counter([trait_avg[0,0], trait_avg[-1,0],trait_avg[0,-1],trait_avg[-1,-1]])

  for cc in c1:
    if c1[cc]>=2:
      trait_avg_nodata = cc
      
  # Repeat the same operation on the 2nd band
  c2 = Counter([trait_std[0,0], trait_std[-1,0],trait_std[0,-1],trait_std[-1,-1]])

  for cc in c2:
    if c2[cc]>=2:
      trait_std_nodata = cc

  
  msk1=None
  msk2=None
  
  if trait_avg_nodata !=9999:
    val_nodata = trait_avg_nodata  #trait_avg[0,0]
    msk1 = (trait_avg == val_nodata)
    
  if trait_std_nodata !=9999:
    val_nodata = trait_std_nodata  #trait_std[0,0]
    msk2 = (trait_std == val_nodata)

  
  if msk1 is not None and msk2 is not None:
    total_mask = (msk1) & (msk2) 
    trait_avg[total_mask] = out_nodata_val
    trait_std[total_mask] = out_nodata_val
  else:
    print("No nodata")
    
  src_ds = None

  return np.concatenate((trait_avg[np.newaxis,:,:],trait_std[np.newaxis,:,:])), src_geotransform, src_wkt
  

# Save output GeoTIFF file, convert it to BSQ, and set nodata value
def save_img(outpath, data_array_bsq, out_geotransform, dest_wkt, data_type=gdal.GDT_Float32):
    
    out_drv = gdal.GetDriverByName( 'GTiff' )
    
    n_band, n_row, n_col = data_array_bsq.shape
    
    dest1 = out_drv.Create(outpath, n_col, \
            n_row, n_band, data_type, options=["INTERLEAVE=BAND","TILED=YES","COMPRESS=LZW"])

    new_geo = out_geotransform
    
    # Set the geotransform
    dest1.SetGeoTransform( new_geo )
    dest1.SetProjection ( dest_wkt )

    for iband in range(n_band):
      dest1.GetRasterBand(iband+1).WriteArray(data_array_bsq[iband,:,:])
      dest1.GetRasterBand(iband+1).SetNoDataValue(out_nodata_val)

    dest1 = None


# main loop for all tif, new tif file has the name "*_v2.tif"    
for trait_file in glob.glob(pat):
  
  print(trait_file)
  
  new_bsq_img , src_geotransform, src_wkt = read_img(trait_file)
  #print(new_bsq_img.shape)
  
  #new_name = trait_file[:-4]+out_suffix
  new_name = os.path.join(outpath,os.path.basename(trait_file)[:-4]+out_suffix)
  
  save_img(new_name, new_bsq_img, src_geotransform, src_wkt)
