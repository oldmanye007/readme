


import argparse,warnings,copy
import numpy as np, os,json, sys

sys.path.append('./HyTools-sandbox')   # need to modify the path


import hytools as ht
from hytools.brdf import *
from hytools.topo_correction import *
from hytools.helpers import *
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
              print(wave_band, abs(wavelengths[wave_band]-single_wave),single_wave,wavelengths[wave_band])
              return {'flag':False}            
            match_index_list = match_index_list + [wave_band]
    return {'flag':True, 'index':match_index_list}
        
    
    
def main():
    '''
    Generate topographic and BRDF correction coefficients. Corrections can be calculated on individual images
    or groups of images.
    '''
    parser = argparse.ArgumentParser(description = "In memory trait mapping tool.")
    parser.add_argument("--img", help="Input image/directory pathname",required=True,nargs = '*', type = str)
    parser.add_argument("--obs", help="Input observables pathname", required=False,nargs = '*', type = str)
    parser.add_argument("--od", help="Ouput directory", required=True, type = str)
    parser.add_argument("--pref", help="Coefficient filename prefix", required=True, type = str)
    parser.add_argument("--brdf", help="Perform BRDF correction",action='store_true')
    parser.add_argument("--kernels", help="Li and Ross kernel types",nargs = 2, type =str)
    parser.add_argument("--topo", help="Perform topographic correction", action='store_true')
    parser.add_argument("--mask", help="Image mask type to use", action='store_true')
    parser.add_argument("--mask_threshold", help="Mask threshold value", nargs = '*', type = float)
    parser.add_argument("--samp_perc", help="Percent of unmasked pixels to sample", type = float,default=1.0)
    parser.add_argument("--agmask", help="ag / urban mask file", required=False, type = str)
    parser.add_argument("--topo_sep", help="In multiple image mode, perform topographic correction in a image-based fasion", action='store_true')
    
    args = parser.parse_args()
    

    if not args.od.endswith("/"):
        args.od+="/"

    if len(args.img) == 1:
        image = args.img[0]
    
        #Load data objects memory
        if image.endswith(".h5"):
            hyObj = ht.openHDF(image,load_obs = True)
        else:
            hyObj = ht.openENVI(image)
            hyObj.load_obs(args.obs[0])
        hyObj.create_bad_bands([[300,400],[1330,1430],[1800,1960],[2450,2600]])
        
        # no data  / ignored values varies by product
        hyObj.no_data =-0.9999
        
        hyObj.load_data()
        
        # Generate mask
        if args.mask:
            ir = hyObj.get_wave(850)
            red = hyObj.get_wave(665)
            ndvi = (1.0*ir-red)/(1.0*ir+red)

            ag_mask=0
            if args.agmask:
              ag_mask =  np.fromfile(args.agmask, dtype=np.uint8).reshape((hyObj.lines,hyObj.columns))

            hyObj.mask = (ndvi > 0.01) & (ndvi < 1.0) & (ir != hyObj.no_data) & (ag_mask ==0)
            del ir,red #,ndvi
        else:
            hyObj.mask = np.ones((hyObj.lines,hyObj.columns)).astype(bool)
            print("Warning no mask specified, results may be unreliable!")

        # Generate cosine i and c1 image for topographic correction
        
        if args.topo:    
            
            topo_coeffs = {}
            topo_coeffs['wavelengths'] = hyObj.wavelengths[hyObj.bad_bands].tolist() 
            topo_coeffs['c'] = []
            cos_i = calc_cosine_i(hyObj.solar_zn, hyObj.solar_az, hyObj.azimuth , hyObj.slope)
            c1 = np.cos(hyObj.solar_zn) 
            c2 = np.cos(hyObj.slope)
            
            terrain_msk = (cos_i > 0.12)  & (hyObj.slope > 0.087)
            topomask = hyObj.mask #& (cos_i > 0.12)  & (hyObj.slope > 0.087) 
              
        # Gernerate scattering kernel images for brdf correction
        if args.brdf:
        
            if args.mask_threshold:
              ndvi_thres = [0.005]+ args.mask_threshold +[1.0]
              total_bin = len(args.mask_threshold)+1
            else:
              ndvi_thres = [0.005,1.0]
              total_bin = 1  

            brdfmask = np.ones(( total_bin, hyObj.lines, hyObj.columns )).astype(bool)
            
            for ibin in range(total_bin):
              brdfmask[ibin,:,:] = hyObj.mask & (ndvi > ndvi_thres[ibin]) & (ndvi <= ndvi_thres[ibin+1]) &  (hyObj.sensor_zn > np.radians(2))

        
            li,ross =  args.kernels
            # Initialize BRDF dictionary
            
            brdf_coeffs_List = [] #initialize
            brdf_mask_stat = np.zeros(total_bin)
            
            for ibin in range(total_bin):
                brdf_mask_stat[ibin] = np.count_nonzero(brdfmask[ibin,:,:])
                if brdf_mask_stat[ibin] < 100:
                  continue
                  
                brdf_coeffs = {}
                brdf_coeffs['li'] = li
                brdf_coeffs['ross'] = ross
                brdf_coeffs['ndvi_lower_bound'] = ndvi_thres[ibin]
                brdf_coeffs['ndvi_upper_bound'] = ndvi_thres[ibin+1]
                brdf_coeffs['wavelengths'] = hyObj.wavelengths[hyObj.bad_bands].tolist() 
                brdf_coeffs['fVol'] = []
                brdf_coeffs['fGeo'] = []
                brdf_coeffs['fIso'] = []
                brdf_coeffs_List.append(brdf_coeffs)
                
            k_vol = generate_volume_kernel(hyObj.solar_az,hyObj.solar_zn,hyObj.sensor_az,hyObj.sensor_zn, ross = ross)
            k_geom = generate_geom_kernel(hyObj.solar_az,hyObj.solar_zn,hyObj.sensor_az,hyObj.sensor_zn,li = li)
            k_finite = np.isfinite(k_vol) & np.isfinite(k_geom)
    
        # Cycle through the bands and calculate the topographic and BRDF correction coefficients
        print("Calculating image correction coefficients.....")
        iterator = hyObj.iterate(by = 'band')
        
        if args.topo or args.brdf:
            while not iterator.complete:   
                band = iterator.read_next() 
                #mask_finite = band > 0.0001
                band_msk = (band>0.001) & (band<0.9)
                progbar(iterator.current_band+1, len(hyObj.wavelengths), 100)
                #Skip bad bands
                if hyObj.bad_bands[iterator.current_band]:
                    # Generate topo correction coefficients
                    if args.topo:
 
                        topomask_b = topomask & band_msk
                            
                        topo_coeff= generate_topo_coeff_band(band,topomask_b & terrain_msk,cos_i)
                        topo_coeffs['c'].append(topo_coeff)

                    # Gernerate BRDF correction coefficients
                    if args.brdf:
                        if args.topo:
                            # Apply topo correction to current band                            
                            correctionFactor = (c2*c1 + topo_coeff  )/(cos_i + topo_coeff)
                            correctionFactor = correctionFactor*topomask_b + 1.0*(1-topomask_b) # only apply to orographic area
                            band = band* correctionFactor
                            
                        for ibin in range(total_bin):
                            if brdf_mask_stat[ibin]<100:
                              continue
                              
                            band_msk_new = (band>0.001) & (band<0.9)  
                              
                            fVol,fGeo,fIso =  generate_brdf_coeff_band(band,brdfmask[ibin,:,:] & band_msk & k_finite & band_msk_new ,k_vol,k_geom)
                            brdf_coeffs_List[ibin]['fVol'].append(fVol)
                            brdf_coeffs_List[ibin]['fGeo'].append(fGeo)
                            brdf_coeffs_List[ibin]['fIso'].append(fIso)
            #print()

           
    # Compute topographic and BRDF coefficients using data from multiple scenes
    elif len(args.img) > 1:
    
        if args.brdf:
        
            li,ross =  args.kernels
            if args.mask_threshold:
              ndvi_thres = [0.005]+ args.mask_threshold +[1.0]
              total_bin = len(args.mask_threshold)+1
            else:
              ndvi_thres = [0.005,1.0]
              total_bin = 1   

            brdf_coeffs_List = [] #initialize
            brdf_mask_stat = np.zeros(total_bin)

            for ibin in range(total_bin):
                brdf_coeffs = {}
                brdf_coeffs['li'] = li
                brdf_coeffs['ross'] = ross
                brdf_coeffs['ndvi_lower_bound'] = ndvi_thres[ibin]
                brdf_coeffs['ndvi_upper_bound'] = ndvi_thres[ibin+1]
                #brdf_coeffs['wavelengths'] = hyObj.wavelengths[hyObj.bad_bands].tolist() 
                brdf_coeffs['fVol'] = []
                brdf_coeffs['fGeo'] = []
                brdf_coeffs['fIso'] = []
                brdf_coeffs_List.append(brdf_coeffs)

            
        hyObj_dict = {}
        sample_dict = {}
        sample_k_vol = []
        sample_k_geom = []
        sample_cos_i = []
        sample_c1 = []
        sample_slope = []
        sample_ndvi = []
        sample_index = [0]
        sub_total_sample_size = 0
        ndvi_mask_dict  = {}
        
        for i,image in enumerate(args.img):
            #Load data objects memory
            if image.endswith(".h5"):
                hyObj = ht.openHDF(image,load_obs = True)
            else:
                hyObj = ht.openENVI(image)
                hyObj.load_obs(args.obs[i])
            hyObj.create_bad_bands([[300,400],[1330,1430],[1800,1960],[2450,2600]])
            hyObj.no_data = -0.9999
            hyObj.load_data()
            
            # Generate mask
            if args.mask:
                ir = hyObj.get_wave(850)
                red = hyObj.get_wave(665)
                ndvi = (1.0*ir-red)/(1.0*ir+red)
                hyObj.mask  =  (ndvi > 0.01) & (ndvi <= 1.0) & (ir != hyObj.no_data)#(ndvi > .5) & (ir != hyObj.no_data) #(ndvi > .01) & (ir != hyObj.no_data)

                del ir,red  #,ndvi
            else:
                hyObj.mask = np.ones((hyObj.lines,hyObj.columns)).astype(bool)
                print("Warning no mask specified, results may be unreliable!")
    
            # Generate sampling mask
            sampleArray = np.zeros(hyObj.mask.shape).astype(bool)
            idx = np.array(np.where(hyObj.mask == True)).T
            
            #np.random.seed(0)  # just for test
            
            idxRand= idx[np.random.choice(range(len(idx)),int(len(idx)*args.samp_perc), replace = False)].T   # actually used
            sampleArray[idxRand[0],idxRand[1]] = True
            sample_dict[i] = sampleArray
            
            print(idxRand.shape)
            sub_total_sample_size+=idxRand.shape[1]
            sample_index = sample_index+ [sub_total_sample_size]
            
            # Initialize and store band iterator
            hyObj_dict[i] = copy.copy(hyObj).iterate(by = 'band')
            
            # Generate cosine i and c1 image for topographic correction
            if args.topo:   
                # Initialize topographic correction dictionary
                topo_coeffs = {}
                topo_coeffs['wavelengths'] = hyObj.wavelengths[hyObj.bad_bands].tolist() 
                topo_coeffs['c'] = []
                sample_cos_i += calc_cosine_i(hyObj.solar_zn, hyObj.solar_az, hyObj.azimuth , hyObj.slope)[sampleArray].tolist()
                sample_c1 += (np.cos(hyObj.solar_zn) * np.cos( hyObj.slope))[sampleArray].tolist()
                sample_slope += (hyObj.slope)[sampleArray].tolist()
            # Gernerate scattering kernel images for brdf correction
            if args.brdf:
            
                sample_ndvi += (ndvi)[sampleArray].tolist()  
                               
                for ibin in range(total_bin):
                    brdf_coeffs_List[ibin]['wavelengths'] = hyObj.wavelengths[hyObj.bad_bands].tolist()     
                
                sample_k_vol += generate_volume_kernel(hyObj.solar_az,hyObj.solar_zn,hyObj.sensor_az,hyObj.sensor_zn, ross = ross)[sampleArray].tolist()
                sample_k_geom += generate_geom_kernel(hyObj.solar_az,hyObj.solar_zn,hyObj.sensor_az,hyObj.sensor_zn,li = li)[sampleArray].tolist()
        
            #del ndvi, topomask, brdfmask
            if args.mask:
                del ndvi
            
        sample_k_vol = np.array(sample_k_vol)
        sample_k_geom = np.array(sample_k_geom)
        sample_cos_i = np.array(sample_cos_i)
        sample_c1= np.array(sample_c1)
        sample_slope = np.array(sample_slope)
        sample_ndvi = np.array(sample_ndvi)
        
        if args.topo_sep:
            topo_coeff_list=[]
            for _ in range(len(args.img)):
                topo_coeff_list+= [{'c':[], 'wavelengths':hyObj.wavelengths[hyObj.bad_bands].tolist()}] 

        if args.brdf:
            for ibin in range(total_bin):
                ndvi_mask = (sample_ndvi>brdf_coeffs_List[ibin]['ndvi_lower_bound']) & (sample_ndvi <=brdf_coeffs_List[ibin]['ndvi_upper_bound'])
                ndvi_mask_dict[ibin] = ndvi_mask
                count_ndvi= np.count_nonzero(ndvi_mask)
                brdf_mask_stat[ibin] = count_ndvi
        
        # Calculate bandwise correction coefficients
        print("Calculating image correction coefficients.....")
        current_progress = 0
        for w,wave in enumerate(hyObj.wavelengths):
            #progbar(current_progress, len(hyObj.wavelengths) * len(args.img), 100)
            wave_samples = []
            for i,image in enumerate(args.img):
                wave_samples +=  hyObj_dict[i].read_next()[sample_dict[i]].tolist()
                current_progress+=1
            
            if hyObj.bad_bands[hyObj_dict[i].current_band]:
                wave_samples = np.array(wave_samples)
                
                # Generate cosine i and c1 image for topographic correction
                if args.topo:    
                  if args.topo_sep==False:
                    #topo_coeff  = generate_topo_coeff_band(wave_samples,(wave_samples>10) & (wave_samples<9000) & (sample_cos_i> 0.12) &  (sample_slope> 0.05) ,sample_cos_i)
                    topo_coeff  = generate_topo_coeff_band(wave_samples,(wave_samples>0.01) & (wave_samples<0.9) & (sample_cos_i> 0.12) &  (sample_slope> 0.087) ,sample_cos_i)
                    topo_coeffs['c'].append(topo_coeff)
                    correctionFactor = (sample_c1 + topo_coeff)/(sample_cos_i + topo_coeff)
                    #print(correctionFactor)
                    wave_samples = wave_samples* correctionFactor
                  else:                   

                    for i in range(len(args.img)):
                        
                        wave_samples_sub = wave_samples[sample_index[i]:sample_index[i+1]]
                        sample_cos_i_sub = sample_cos_i[sample_index[i]:sample_index[i+1]]
                        sample_slope_sub = sample_slope[sample_index[i]:sample_index[i+1]]
                        sample_c1_sub = sample_c1[sample_index[i]:sample_index[i+1]]

                        
                        #topo_coeff  = generate_topo_coeff_band(wave_samples,(wave_samples_sub>10) & (wave_samples_sub<9000) & (sample_cos_i_sub> 0.12) &  (sample_slope_sub> 0.05) ,sample_cos_i_sub)
                        topo_coeff  = generate_topo_coeff_band(wave_samples_sub,(wave_samples_sub>0.01) & (wave_samples_sub<0.9) & (sample_cos_i_sub> 0.12) &  (sample_slope_sub> 0.087) ,sample_cos_i_sub)
                        topo_coeff_list[i]['c'].append(topo_coeff)

                        #topo_coeff_list+=[topo_coeffs_img]
                        correctionFactor = (sample_c1_sub + topo_coeff)/(sample_cos_i_sub + topo_coeff)
                        #print(correctionFactor)
                        wave_samples[sample_index[i]:sample_index[i+1]] = wave_samples_sub* correctionFactor                    
                # Gernerate scattering kernel images for brdf correction
                if args.brdf:
                
                    temp_mask = (wave_samples>0.01) & (wave_samples<0.9) & np.isfinite(sample_k_vol) & np.isfinite(sample_k_geom)
                    temp_mask = temp_mask & (sample_cos_i> 0.12) &  (sample_slope> 0.03) 
                
                    for ibin in range(total_bin):

                        if brdf_mask_stat[ibin]<100:
                          continue
                          
                        fVol,fGeo,fIso = generate_brdf_coeff_band(wave_samples, temp_mask & ndvi_mask_dict[ibin]  ,sample_k_vol,sample_k_geom)    
                        brdf_coeffs_List[ibin]['fVol'].append(fVol)
                        brdf_coeffs_List[ibin]['fGeo'].append(fGeo)
                        brdf_coeffs_List[ibin]['fIso'].append(fIso)
                

    
    
    # Export coefficients to JSON           
    if args.topo:
      if args.topo_sep==False:
        topo_json = "%s%s_topo_coeffs.json" % (args.od,args.pref)
        with open(topo_json, 'w') as outfile:
            json.dump(topo_coeffs,outfile) 
      else:
        #print(len(topo_coeff_list))
        for i_img in range(len(args.img)):
            topo_json = "%s%s_%s_topo_coeffs.json" % (args.od,args.pref,i_img)
            with open(topo_json, 'w') as outfile:
                json.dump(topo_coeff_list[i_img],outfile)            


    if args.brdf:
      if total_bin>0:
        for ibin in range(total_bin):            
            if brdf_mask_stat[ibin]<100:
              continue
            brdf_json = "%s%s_brdf_coeffs_%s.json" % (args.od,args.pref,str(ibin+1))
            with open(brdf_json, 'w') as outfile:
                json.dump(brdf_coeffs_List[ibin],outfile)
      else:
            brdf_json = "%s%s_brdf_coeffs_1.json" % (args.od,args.pref)
            with open(brdf_json, 'w') as outfile:
                json.dump(brdf_coeffs,outfile)                

                
if __name__== "__main__":
    main()

#python topo_brdf_coeffs.py --img --pref --brdf --kernels dense thick --mask --mask_threshold .7 --samp_perc 0.1


