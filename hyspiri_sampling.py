import argparse,gdal,warnings
import numpy as np, os,pandas as pd,glob, sys

# need to modify the path
sys.path.append('./HyTools-sandbox')

import hytools as ht
from hytools.sampling import *

home = os.path.expanduser("~")

warnings.filterwarnings("ignore")

def main():
    '''
    Perform sampling with text point.
    '''
    
    parser = argparse.ArgumentParser(description = "Sampling / averaging tool.")
    parser.add_argument("-img", help="Input image pathname",required=True, type = str)
    parser.add_argument("-od", help="Output directory for all resulting products", required=True, type = str)
    parser.add_argument("-pnt", help="Point location files with UID", required=True, type = str)
    #parser.add_argument("--obs", help="Input observables pathname", required=False, type = str)

    args = parser.parse_args()
    
    points_file = args.pnt
    print(points_file)

    #Load data objects memory
    if args.img.endswith(".h5"):
        hyObj = ht.openHDF(args.img,load_obs = False)
    else:
        hyObj = ht.openENVI(args.img)

    if not args.od.endswith("/"):
        args.od+="/"
    #hyObj.create_bad_bands([[300,400],[1330,1430],[1800,1960],[2450,2600]])
    #hyObj.create_bad_bands([[300,420],[1320,1430],[1800,1960],[2450,2600]])   
    
    #hyObj.load_obs(args.obs)    
    hyObj.no_data =-0.9999
    
    hyObj.load_data()
 
    #out_df = point.point2spec(hyObj, points_file, 'Plot', 'x_utm11n', 'y_utm11n', 32611, n_neighbor=8, use_band_list=False, band_list=[])  #4326  
    out_df = point.point2spec(hyObj, points_file, 'Plot', 'lon', 'lat', 4326, n_neighbor=8, use_band_list=False, band_list=[])  #4326  
    
    if not out_df is None:
        img_row = out_df['img_row'].values.astype(np.int)
        img_col = out_df['img_col'].values.astype(np.int)
        
        #out_obs_df = obs_point2spec(hyObj, img_row, img_col)
        
        img_base_name=os.path.basename(args.img)
        #out_obs_df.to_csv(args.od+img_base_name+"_obs_df.csv",index=False)

        out_df.to_csv(args.od+img_base_name+"_spec_df.csv",index=False) #
    


if __name__== "__main__":
    main()
