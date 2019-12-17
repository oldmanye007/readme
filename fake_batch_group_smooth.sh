
#echo $1 # rfl image dir
#echo $2 # obs file dir
#echo $3 #file list
#echo $4 # output dir
#echo $5 # coeffs dir 
#echo $6 # group tag 

str_img=''
str_obs=''


# Load images names from file list and merge to two variable lists
while IFS= read -r imgbase
do
    #echo $imgbase

    imgname=$(ls $1/$imgbase*rfl* | head -1)
    obsname=$(ls $2/$imgbase*obs* | head -1)
    
    #echo "$imgname"
    #echo "$obsname"

    # join the variable lists
    str_img="$str_img $imgname"
    str_obs="$str_obs $obsname"

    
done < $3

# Estimate BRDF correction coefficients in group mode. The two variable lists are gathered in the previous loop
python topo_brdf_coeffs_modify_hybeta4.py  --img $str_img  --obs $str_obs  --od $4  --pref $6   --topo --brdf  --mask --mask_threshold 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 --kernels sparse thick  --samp_perc 0.1 --topo_sep



# Perform TOPO+BRDF correction image by image and output desired product. Interpolation smoothing for BRDF coefficients is used in this example.
while IFS= read -r imgbase
do
    #echo $imgbase

    imgname=$(ls $1/$imgbase*rfl* | head -1)
    obsname=$(ls $2/$imgbase*obs* | head -1)
    
    #echo "$imgname"
    #echo "$obsname"
    
    python image_to_traits_modify_hybeta4.py -img $imgname --obs $obsname -od $4 --mask --mask_threshold 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 --brdf $4/$6 --topo $4/"$imgbase"_topo_coeffs.json   -coeffs $5 -smooth I

done < $3


#bash fake_batch_group_smooth.sh /hyspiri_image_dir/ /hyspiri_obs_dir/ /hyspiri_file_list_dir/fake_list.txt /hyspiri/coeff_out_dir/ /hyspiri/out_traits_dir/ hyspiri_group_f130626_18bins