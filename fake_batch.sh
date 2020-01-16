
#echo $1 # rfl image dir
#echo $2 # obs file dir
#echo $3 #file list
#echo $4 # output dir
#echo $5 # coeffs dir 

while IFS= read -r imgbase
do
    #echo $imgbase

    imgname=$(ls $1/$imgbase*rfl* | head -1)
    obsname=$(ls $2/$imgbase*obs* | head -1)
    
    echo "$imgname"
    echo "$obsname"
    
    #python topo_brdf_coeffs_modify_hybeta4.py  --img $imgname  --obs $obsname  --od $4  --pref $imgbase    --kernels sparse thick  --mask  --mask_threshold 0.3 0.7  --brdf --topo
    #python image_to_traits_modify_hybeta4.py -img $imgname --obs $obsname -od $4 --mask --mask_threshold 0.3 0.7 --brdf $4/$imgbase --topo $4/"$imgbase"_topo_coeffs.json   -coeffs $5
    
done < $3
