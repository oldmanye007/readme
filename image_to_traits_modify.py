import argparse,gdal,copy,sys,warnings
import numpy as np, os,pandas as pd,glob,json, sys

sys.path.append('./HyTools-sandbox')   # need to modify the path

import hytools as ht
from hytools.brdf import *
from hytools.topo_correction import *
from hytools.helpers import *
from hytools.preprocess.resampling import *
from hytools.preprocess.vector_norm import *
from hytools.file_io import array_to_geotiff,writeENVI
home = os.path.expanduser("~")

warnings.filterwarnings("ignore")

def progbar(curr, total, full_progbar):
    frac = curr/total
    filled_progbar = round(frac*full_progbar)
    print('\r', '#'*filled_progbar + '-'*(full_progbar-filled_progbar), '[{:>7.2%}]'.format(frac), end='')

def check_wave_match(hyObj, dst_wave_list):
    
    wavelengths = hyObj.wavelengths[hyObj.bad_bands]
    match_index_list = []
    for single_wave in dst_wave_list:
    
        if hyObj.wavelength_units == "micrometers" and single_wave > 3:
            single_wave/= 1000
        if hyObj.wavelength_units == "nanometers" and single_wave < 3:
            single_wave*= 1000    
    
        if single_wave in wavelengths:
            wave_band = np.argwhere(single_wave == wavelengths)[0][0]
        elif (single_wave  > wavelengths.max()) | (single_wave  < wavelengths.min()):
            return {'flag':False}
        else: 
            wave_band = np.argmin(np.abs(wavelengths - single_wave))    
            if abs(wavelengths[wave_band]-single_wave) > 0.5:
              #print(wave_band, abs(wavelengths[wave_band]-single_wave),single_wave,wavelengths[wave_band])
              return {'flag':False}            
            match_index_list = match_index_list + [wave_band]
    return {'flag':True, 'index':match_index_list}        
    
def main():
    '''
    Perform in-memory trait estimation.
    '''
    
    parser = argparse.ArgumentParser(description = "In memory trait mapping tool.")
    parser.add_argument("-img", help="Input image pathname",required=True, type = str)
    parser.add_argument("--obs", help="Input observables pathname", required=False, type = str)
    parser.add_argument("--out", help="Output full corrected image", required=False, type = str)
    parser.add_argument("-od", help="Output directory for all resulting products", required=True, type = str)
    parser.add_argument("--brdf", help="Perform BRDF correction",type = str, default = '')
    parser.add_argument("--topo", help="Perform topographic correction", type = str, default = '')
    parser.add_argument("--mask", help="Image mask type to use", action='store_true')
    parser.add_argument("--mask_threshold", help="Mask threshold value", nargs = '*', type = float)
    parser.add_argument("--rgbim", help="Export RGBI +Mask image.", action='store_true')
    #parser.add_argument("-coeffs", help="Trait coefficients directory", required=True, type = str)
    parser.add_argument("-coeffs", help="Trait coefficients directory", required=False, type = str)
    parser.add_argument("-nodata", help="New value to assign for no_data values", required=False, type = float, default =-9999)
    args = parser.parse_args()


    traits = glob.glob("%s/*.json" % args.coeffs)
    
    #Load data objects memory
    if args.img.endswith(".h5"):
        hyObj = ht.openHDF(args.img,load_obs = True)
    else:
        hyObj = ht.openENVI(args.img)
    if (len(args.topo) != 0) | (len(args.brdf) != 0):
        hyObj.load_obs(args.obs)
    if not args.od.endswith("/"):
        args.od+="/"
    hyObj.create_bad_bands([[300,400],[1330,1430],[1800,1960],[2450,2600]])
    
    # no data  / ignored values varies by product
    hyObj.no_data =-0.9999
    
    hyObj.load_data()
    
    # Generate mask
    if args.mask:
        ir = hyObj.get_wave(850)
        red = hyObj.get_wave(665)
        ndvi = (ir-red)/(ir+red)
        #mask = (ndvi > args.mask_threshold) & (ir != hyObj.no_data)
        #mask = (ndvi > 0.05) & (ir != hyObj.no_data)
        hyObj.mask = (ndvi > 0.05) & (ir != hyObj.no_data) 
        del ir,red #,ndvi
    else:
        hyObj.mask = np.ones((hyObj.lines,hyObj.columns)).astype(bool)
        print("Warning no mask specified, results may be unreliable!")

    # Generate cosine i and c1 image for topographic correction
    if len(args.topo) != 0:
        with open( args.topo) as json_file:  
            topo_coeffs = json.load(json_file)
            
        topo_coeffs['c'] = np.array(topo_coeffs['c'])   
        cos_i =  calc_cosine_i(hyObj.solar_zn, hyObj.solar_az, hyObj.azimuth , hyObj.slope)
        c1 = np.cos(hyObj.solar_zn)
        c2 = np.cos(hyObj.slope)
        
        topomask = hyObj.mask & (cos_i > 0.12)  & (hyObj.slope > 0.087)
           
    #total_bin = len(args.mask_threshold)+1
    #brdf_coeffs_List = []
    
    # Gernerate scattering kernel images for brdf correction
    if len(args.brdf) != 0:
        brdf_coeffs_List = []

        if (args.mask_threshold):
          total_bin = len(args.mask_threshold)+1
          ndvi_thres = [0.05]+ args.mask_threshold +[1.0]
        else:
          total_bin=1
          ndvi_thres = [0.05, 1.0]
    
        #ndvi_thres = [0.05]+ args.mask_threshold +[1.0]
        brdfmask = np.ones((total_bin, hyObj.lines,hyObj.columns )).astype(bool)
            
        for ibin in range(total_bin):
        
          with open(args.brdf+'_brdf_coeffs_'+str(ibin+1)+'.json') as json_file:  
            brdf_coeffs = json.load(json_file)
            brdf_coeffs['fVol'] = np.array(brdf_coeffs['fVol'])
            brdf_coeffs['fGeo'] = np.array(brdf_coeffs['fGeo'])
            brdf_coeffs['fIso'] = np.array(brdf_coeffs['fIso'])
            brdf_coeffs_List.append(brdf_coeffs)
            
          brdfmask[ibin, :,:] = hyObj.mask & (ndvi > ndvi_thres[ibin]) & (ndvi <= ndvi_thres[ibin+1])            

        
        k_vol = generate_volume_kernel(hyObj.solar_az,hyObj.solar_zn,hyObj.sensor_az,hyObj.sensor_zn, ross = brdf_coeffs_List[0]['ross'])
        k_geom = generate_geom_kernel(hyObj.solar_az,hyObj.solar_zn,hyObj.sensor_az,hyObj.sensor_zn,li = brdf_coeffs_List[0]['li'])
        k_vol_nadir = generate_volume_kernel(hyObj.solar_az,hyObj.solar_zn,hyObj.sensor_az,0, ross = brdf_coeffs_List[0]['ross'])
        k_geom_nadir = generate_geom_kernel(hyObj.solar_az,hyObj.solar_zn,hyObj.sensor_az,0,li = brdf_coeffs_List[0]['li'])
        
        


    if len(traits)!=0:
      
        #Cycle through the chunks and apply topo, brdf, vnorm,resampling and trait estimation steps
        print("Calculating values for %s traits....." % len(traits))
        
        
        # Cycle through trait models and gernerate resampler
        trait_waves_all = []
        trait_fwhm_all = []
        
        for i,trait in enumerate(traits):
            with open(trait) as json_file:  
                trait_model = json.load(json_file)
             
            # Check if wavelength units match
            if trait_model['wavelength_units'] == 'micrometers':
                trait_wave_scaler = 10**3
            else:
                trait_wave_scaler = 1    
            
            # Get list of wavelengths to compare against image wavelengths
            if len(trait_model['vector_norm_wavelengths']) == 0:
                trait_waves_all += list(np.array(trait_model['model_wavelengths'])*trait_wave_scaler)
            else:
                trait_waves_all += list(np.array(trait_model['vector_norm_wavelengths'])*trait_wave_scaler)
            
            trait_fwhm_all += list(np.array(trait_model['fwhm'])* trait_wave_scaler)        
               
        # List of all unique pairs of wavelengths and fwhm    
        trait_waves_fwhm = list(set([x for x in zip(trait_waves_all,trait_fwhm_all)]))
        trait_waves_fwhm.sort(key = lambda x: x[0])
        
        # Create a single set of resampling coefficients for all wavelength and fwhm combos
        #resampling_coeffs = est_transform_matrix(hyObj.wavelengths[hyObj.bad_bands],[x for (x,y) in trait_waves_fwhm] ,hyObj.fwhm[hyObj.bad_bands],[y for (x,y) in trait_waves_fwhm],1)

        # if wavelengths match, no need to resample
        check_wave_match_result = check_wave_match(hyObj, [x for (x,y) in trait_waves_fwhm])
        if (check_wave_match_result['flag']):
            match_flag = True
        else:
            match_flag = False
            resampling_coeffs = est_transform_matrix(hyObj.wavelengths[hyObj.bad_bands],[x for (x,y) in trait_waves_fwhm] ,hyObj.fwhm[hyObj.bad_bands],[y for (x,y) in trait_waves_fwhm],1) #2        
        
    #else:
    #  print('no trait')
      
    hyObj.wavelengths = hyObj.wavelengths[hyObj.bad_bands]
    
    pixels_processed = 0
    iterator = hyObj.iterate(by = 'chunk',chunk_size = (32,hyObj.columns))

    while not iterator.complete:  
        chunk = iterator.read_next()  
        chunk_nodata_mask = chunk[:,:,10] == hyObj.no_data  # 10th band
        pixels_processed += chunk.shape[0]*chunk.shape[1]
        progbar(pixels_processed, hyObj.columns*hyObj.lines, 100)

        # Chunk Array indices
        line_start =iterator.current_line 
        line_end = iterator.current_line + chunk.shape[0]
        col_start = iterator.current_column
        col_end = iterator.current_column + chunk.shape[1]
        
        # Apply TOPO correction 
        if len(args.topo) != 0:
            cos_i_chunk = cos_i[line_start:line_end,col_start:col_end]
            c1_chunk = c1[line_start:line_end,col_start:col_end]
            c2_chunk = c2[line_start:line_end,col_start:col_end]
            topomask_chunk = topomask[line_start:line_end,col_start:col_end,np.newaxis]
            correctionFactor = (c2_chunk[:,:,np.newaxis]*c1_chunk[:,:,np.newaxis] +topo_coeffs['c']  )/(cos_i_chunk[:,:,np.newaxis] + topo_coeffs['c'])
            correctionFactor = correctionFactor*topomask_chunk + 1.0*(1-topomask_chunk)
            chunk = chunk[:,:,hyObj.bad_bands]* correctionFactor
        else:
            chunk = chunk[:,:,hyObj.bad_bands] *1
        
        # Apply BRDF correction 
        if len(args.brdf) != 0:
            # Get scattering kernel for chunks
            k_vol_nadir_chunk = k_vol_nadir[line_start:line_end,col_start:col_end]
            k_geom_nadir_chunk = k_geom_nadir[line_start:line_end,col_start:col_end]
            k_vol_chunk = k_vol[line_start:line_end,col_start:col_end]
            k_geom_chunk = k_geom[line_start:line_end,col_start:col_end]
            
            #veg_mask = brdfmask[:,line_start:line_end,col_start:col_end]
        
            n_wavelength = brdf_coeffs_List[0]['fVol'].shape[0]
            new_k_vol = np.zeros((chunk.shape[0],chunk.shape[1],n_wavelength))
            new_k_geom = np.zeros((chunk.shape[0],chunk.shape[1],n_wavelength))
            new_k_iso = np.zeros((chunk.shape[0],chunk.shape[1],n_wavelength))
            
            #v_msk1 = brdfmask[0,line_start:line_end,col_start:col_end][:,:,np.newaxis]   #| (veg_total == 0)
            #v_msk2 = brdfmask[1,line_start:line_end,col_start:col_end][:,:,np.newaxis]  
            #v_msk3 = brdfmask[2,line_start:line_end,col_start:col_end][:,:,np.newaxis]  
            
            for ibin in range(total_bin):
              
              veg_mask = brdfmask[ibin,line_start:line_end,col_start:col_end][:,:,np.newaxis]
              
              new_k_vol +=  brdf_coeffs_List[ibin]['fVol'] * veg_mask
              new_k_geom += brdf_coeffs_List[ibin]['fGeo'] * veg_mask
              new_k_iso += brdf_coeffs_List[ibin]['fIso'] * veg_mask
            
            #new_k_vol = brdf_df1.k_vol.values * v_msk1 + brdf_df2.k_vol.values * v_msk2 + brdf_df3.k_vol.values * v_msk3  # + brdf_df4.k_vol.values * v_msk4
        
            #new_k_geom = brdf_df1.k_geom.values * v_msk1 + brdf_df2.k_geom.values * v_msk2 + brdf_df3.k_geom.values * v_msk3    #+ brdf_df4.k_geom.values * v_msk4
        
            #new_k_iso = brdf_df1.k_iso.values * v_msk1 + brdf_df2.k_iso.values * v_msk2 + brdf_df3.k_iso.values * v_msk3
            
    
            # Apply brdf correction 
            # eq 5. Weyermann et al. IEEE-TGARS 2015)
            #brdf = np.einsum('i,jk-> jki', brdf_coeffs['fVol'],k_vol_chunk) + np.einsum('i,jk-> jki', brdf_coeffs['fGeo'],k_geom_chunk)  + brdf_coeffs['fIso']
            #brdf_nadir = np.einsum('i,jk-> jki', brdf_coeffs['fVol'],k_vol_nadir_chunk) + np.einsum('i,jk-> jki', brdf_coeffs['fGeo'],k_geom_nadir_chunk)  +brdf_coeffs['fIso']
            #correctionFactor = brdf_nadir/brdf
            #chunk= chunk* correctionFactor
            
            
            
            brdf = np.einsum('ijk,ij-> ijk', new_k_vol,k_vol_chunk) + np.einsum('ijk,ij-> ijk', new_k_geom,k_geom_chunk)  + new_k_iso
            brdf_nadir = np.einsum('ijk,ij-> ijk', new_k_vol,k_vol_nadir_chunk) + np.einsum('ijk,ij-> ijk', new_k_geom,k_geom_nadir_chunk)  + new_k_iso
            correctionFactor = brdf_nadir/brdf  #*veg_total+(1.0-veg_total)
            correctionFactor[brdf == 0.0] = 1.0
            chunk = chunk   * correctionFactor
            
        
        #Reassign no data values
        chunk[chunk_nodata_mask,:] = args.nodata
        

        if len(traits)>0:
        # Resample chunk 
          #chunk_r = np.dot(chunk, resampling_coeffs) 
          if match_flag==False:            
              chunk_r = np.dot(chunk, resampling_coeffs) 
          # subset of chunk
          else:
              chunk_r = chunk[:,:,check_wave_match_result['index']]
                
        # Export RGBIM image
        if args.rgbim:
            dstFile = args.od + os.path.splitext(os.path.basename(args.img))[0] + '_rgbim.tif'
            if line_start + col_start == 0:
                driver = gdal.GetDriverByName("GTIFF")
                #tiff = driver.Create(dstFile,hyObj.columns,hyObj.lines,5,gdal.GDT_Float32)
                tiff = driver.Create(dstFile,hyObj.columns,hyObj.lines,8,gdal.GDT_Float32)
                tiff.SetGeoTransform(hyObj.transform)
                tiff.SetProjection(hyObj.projection)
                for band in range(1,9):
                    tiff.GetRasterBand(band).SetNoDataValue(0)
                tiff.GetRasterBand(8).WriteArray(hyObj.mask)

                del tiff,driver
            # Write rgbi chunk
            rgbi_geotiff = gdal.Open(dstFile, gdal.GA_Update)
            
            for i,wave in enumerate([480,560,660,850,976,1650,2217],start=1):
                    #print(wave, band)
                    band = hyObj.wave_to_band(wave)
                    
                    rgbi_geotiff.GetRasterBand(i).WriteArray(chunk[:,:,band], col_start, line_start)
            rgbi_geotiff = None
        
        # Export BRDF and topo corrected image
        if args.out:
            if line_start + col_start == 0:
                #output_name = args.od + os.path.splitext(os.path.basename(args.img))[0] + "_topo_brdf" 
                output_name = args.od + os.path.splitext(os.path.basename(args.img))[0] + args.out 
                header_dict =hyObj.header_dict
                # Update header
                header_dict['wavelength']= header_dict['wavelength'][hyObj.bad_bands]
                header_dict['fwhm'] = header_dict['fwhm'][hyObj.bad_bands]
                #header_dict['bbl'] = header_dict['bbl'][hyObj.bad_bands]
                if 'band names' in header_dict:  
                  del header_dict['band names']
                header_dict['bands'] = int(hyObj.bad_bands.sum())
                
                # clean ENVI header
                header_dict.pop('band names', None)
                header_dict.pop('correction factors', None)
                header_dict.pop('bbl', None)                
                header_dict.pop('smoothing factors', None)               
                
                writer = writeENVI(output_name,header_dict)
            writer.write_chunk(chunk,iterator.current_line,iterator.current_column)
            if iterator.complete:
                writer.close()
                
        for i,trait in enumerate(traits):
            dstFile = args.od + os.path.splitext(os.path.basename(args.img))[0] +"_" +os.path.splitext(os.path.basename(trait))[0] +".tif"
            
            # Trait estimation preparation
            if line_start + col_start == 0:
                
                with open(trait) as json_file:  
                    trait_model = json.load(json_file)
       
                intercept = np.array(trait_model['intercept'])
                coefficients = np.array(trait_model['coefficients'])
                transform = trait_model['transform']
                
                # Get list of wavelengths to compare against image wavelengths
                if len(trait_model['vector_norm_wavelengths']) == 0:
                    dst_waves = np.array(trait_model['model_wavelengths'])*trait_wave_scaler
                else:
                    dst_waves = np.array(trait_model['vector_norm_wavelengths'])*trait_wave_scaler
                
                dst_fwhm = np.array(trait_model['fwhm'])* trait_wave_scaler
                model_waves = np.array(trait_model['model_wavelengths'])* trait_wave_scaler
                model_fwhm = [dict(zip(dst_waves, dst_fwhm))[x] for x in model_waves]
                
                vnorm_band_mask = [x in zip(dst_waves,dst_fwhm) for x in trait_waves_fwhm]
                model_band_mask = [x in zip(model_waves,model_fwhm) for x in trait_waves_fwhm]
                
                vnorm_band_mask = np.array(vnorm_band_mask)  # convert list to numpy array, otherwise True/False will be treated as 1/0, which is the 2nd/1st band 
                model_band_mask = np.array(model_band_mask)   # convert list to numpy array, otherwise True/False will be treated as 1/0, which is the 2nd/1st band

                if trait_model['vector_norm']:
                    vnorm_scaler = trait_model["vector_scaler"]
                else:
                    vnorm_scaler = None

                # Initialize trait dictionary
                if i == 0:
                    trait_dict = {}
                trait_dict[i] = [coefficients,intercept,trait_model['vector_norm'],vnorm_scaler,vnorm_band_mask,model_band_mask,transform]
        
                # Create geotiff driver
                driver = gdal.GetDriverByName("GTIFF")
                tiff = driver.Create(dstFile,hyObj.columns,hyObj.lines,2,gdal.GDT_Float32)
                tiff.SetGeoTransform(hyObj.transform)
                tiff.SetProjection(hyObj.projection)
                tiff.GetRasterBand(1).SetNoDataValue(0)
                tiff.GetRasterBand(2).SetNoDataValue(0)
                del tiff,driver
            
            coefficients,intercept,vnorm,vnorm_scaler,vnorm_band_mask,model_band_mask,transform = trait_dict[i]

            chunk_t =np.copy(chunk_r)

            if vnorm:                    
                chunk_t[:,:,vnorm_band_mask] = vector_normalize_chunk(chunk_t[:,:,vnorm_band_mask],vnorm_scaler)

            
            if transform == "log(1/R)":
                chunk_t[:,:,model_band_mask] = np.log(1/chunk_t[:,:,model_band_mask] )
                

            trait_mean,trait_std = apply_plsr_chunk(chunk_t[:,:,model_band_mask],coefficients,intercept)
            
            
            # Change no data pixel values
            trait_mean[chunk_nodata_mask] = args.nodata
            trait_std[chunk_nodata_mask] = args.nodata

            # Write trait estimate to file
            trait_geotiff = gdal.Open(dstFile, gdal.GA_Update)
            trait_geotiff.GetRasterBand(1).WriteArray(trait_mean, col_start, line_start)
            trait_geotiff.GetRasterBand(2).WriteArray(trait_std, col_start, line_start)
            trait_geotiff = None

if __name__== "__main__":
    main()

