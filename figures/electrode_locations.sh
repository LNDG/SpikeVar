# #!/bin/bash
# This script creates the electrode masks for all participants      

# source exploconfig.sh
root="/Users/kamp/PhD/spikevar/repo/SpikeVar"
figure_dir="${root}/figures"
position_file="${figure_dir}/electrodePos.csv"
mni_file="${figure_dir}/MNI152_T1_1mm_brain.nii.gz"
sphere_size=2

fslmaths $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz -mul 0 ./overlays/amygdala_mask # initiating combined mask across all subjects
fslmaths $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz -mul 0 ./overlays/hpc_mask 

while IFS=';' read id l_amygdala r_amygdala l_hpc r_hpc
do    
    if [ "${id:0:1}" != "P" ] || [ "${id}" = "Patient ID" ]; then
        continue
    fi
    echo "Patient $id"
    fslmaths $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz -mul 0 ./overlays/${id}_amygdala_mask
    # left amygdala
    if [[ ${l_amygdala:0:3} != "n/a" ]]; then   
        x_lamygdala=$(echo $l_amygdala | LC_ALL=C cut -d',' -f 1 | LC_ALL=C grep -o -E '[-0-9]+.[0-9]+')  # get x coordinate and remove trailing characters
        y_lamygdala=$(echo $l_amygdala | LC_ALL=C cut -d',' -f 2 | LC_ALL=C grep -o -E '[-0-9]+.[0-9]+')
        z_lamygdala=$(echo $l_amygdala | LC_ALL=C cut -d',' -f 3 | LC_ALL=C grep -o -E '[-0-9]+.[0-9]+')
        echo "left amygdala $x_lamygdala $y_lamygdala $z_lamygdala"        

        x_lamygdala=$(echo "($x_lamygdala * -1 + 90)/2" | LC_ALL=C bc) # convert to voxel coordinates
        y_lamygdala=$(echo "($y_lamygdala * 1 + 126)/2" | LC_ALL=C bc)
        z_lamygdala=$(echo "($z_lamygdala * 1 + 72)/2" | LC_ALL=C bc)

        fslmaths $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz -mul 0 -add 1 -roi $x_lamygdala 1 $y_lamygdala 1 $z_lamygdala 1 0 1 ./overlays/${id}_leftamygdala -odt float
        fslmaths ./overlays/${id}_leftamygdala.nii.gz -kernel sphere $sphere_size -fmean ./overlays/${id}_leftamygdala_sphere -odt float
        fslmaths ./overlays/${id}_leftamygdala_sphere.nii.gz -bin ./overlays/${id}_leftamygdala_bin.nii.gz  
        fslmaths ./overlays/${id}_amygdala_mask -add ./overlays/${id}_leftamygdala_bin.nii.gz ./overlays/${id}_amygdala_mask
        rm ./overlays/${id}_leftamygdala.nii.gz ./overlays/${id}_leftamygdala_sphere.nii.gz ./overlays/${id}_leftamygdala_bin.nii.gz # delete files that we don't need anymore
    fi

    # right amygdala
    if [[ ${r_amygdala:0:3} != "n/a" ]]; then   
        x_ramygdala=$(echo $r_amygdala | LC_ALL=C cut -d ',' -f 1 | LC_ALL=C grep -o -E '[-0-9]+.[0-9]+')  
        y_ramygdala=$(echo $r_amygdala | LC_ALL=C cut -d ',' -f 2 | LC_ALL=C grep -o -E '[-0-9]+.[0-9]+')
        z_ramygdala=$(echo $r_amygdala | LC_ALL=C cut -d ',' -f 3 | LC_ALL=C grep -o -E '[-0-9]+.[0-9]+')
        echo "right amygdala $x_ramygdala $y_ramygdala $z_ramygdala"

        x_ramygdala=$(echo "($x_ramygdala * -1 + 90)/2" | bc) # convert to voxel coordinates
        y_ramygdala=$(echo "($y_ramygdala * 1 + 126)/2" | bc)
        z_ramygdala=$(echo "($z_ramygdala * 1 + 72)/2" | bc)

        fslmaths $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz -mul 0 -add 1 -roi $x_ramygdala 1 $y_ramygdala 1 $z_ramygdala 1 0 1 ./overlays/${id}_rightamygdala -odt float
        fslmaths ./overlays/${id}_rightamygdala.nii.gz -kernel sphere $sphere_size -fmean ./overlays/${id}_rightamygdala_sphere -odt float
        fslmaths ./overlays/${id}_rightamygdala_sphere.nii.gz -bin ./overlays/${id}_rightamygdala_bin.nii.gz  
        fslmaths ./overlays/${id}_amygdala_mask -add ./overlays/${id}_rightamygdala_bin.nii.gz ./overlays/${id}_amygdala_mask
        rm ./overlays/${id}_rightamygdala.nii.gz ./overlays/${id}_rightamygdala_sphere.nii.gz ./overlays/${id}_rightamygdala_bin.nii.gz
    fi

    fslmaths $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz -mul 0 ./overlays/${id}_hpc_mask
    # left hippocampus
    if [[ ${l_hpc:0:3} != "n/a" ]]; then   
        x_lhpc=$(echo $l_hpc | LC_ALL=C cut -d',' -f 1 | LC_ALL=C grep -o -E '[-0-9]+.[0-9]+')  
        y_lhpc=$(echo $l_hpc | LC_ALL=C cut -d',' -f 2 | LC_ALL=C grep -o -E '[-0-9]+.[0-9]+')
        z_lhpc=$(echo $l_hpc | LC_ALL=C cut -d',' -f 3 | LC_ALL=C grep -o -E '[-0-9]+.[0-9]+')
        echo "left hpc $x_lhpc $y_lhpc $z_lhpc"

        x_lhpc=$(echo "($x_lhpc * -1 + 90)/2" | bc) # convert to voxel coordinates; divide by 2 if you're using the MNI 2mm map
        y_lhpc=$(echo "($y_lhpc * 1 + 126)/2" | bc)
        z_lhpc=$(echo "($z_lhpc * 1 + 72)/2" | bc)

        fslmaths $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz -mul 0 -add 1 -roi $x_lhpc 1 $y_lhpc 1 $z_lhpc 1 0 1 ./overlays/${id}_lefthpc -odt float
        fslmaths ./overlays/${id}_lefthpc.nii.gz -kernel sphere $sphere_size -fmean ./overlays/${id}_lefthpc_sphere -odt float
        fslmaths ./overlays/${id}_lefthpc_sphere.nii.gz -bin ./overlays/${id}_lefthpc_bin.nii.gz  
        fslmaths ./overlays/${id}_hpc_mask -add ./overlays/${id}_lefthpc_bin.nii.gz ./overlays/${id}_hpc_mask
        rm ./overlays/${id}_lefthpc.nii.gz ./overlays/${id}_lefthpc_sphere.nii.gz ./overlays/${id}_lefthpc_bin.nii.gz
    fi

    # right hippocampus
    if [[ ${r_hpc:0:3} != "n/a" ]]; then   
        x_rhpc=$(echo $r_hpc | LC_ALL=C cut -d',' -f 1 | LC_ALL=C grep -o -E '[-0-9]+.[0-9]+')  
        y_rhpc=$(echo $r_hpc | LC_ALL=C cut -d',' -f 2 | LC_ALL=C grep -o -E '[-0-9]+.[0-9]+')
        z_rhpc=$(echo $r_hpc | LC_ALL=C cut -d',' -f 3 | LC_ALL=C grep -o -E '[-0-9]+.[0-9]+')
        echo "right hpc $x_rhpc $y_rhpc $z_rhpc"
        
        x_rhpc=$(echo "($x_rhpc * -1 + 90)/2" | LC_ALL=C bc) # convert to voxel coordinates
        y_rhpc=$(echo "($y_rhpc * 1 + 126)/2" | LC_ALL=C bc)
        z_rhpc=$(echo "($z_rhpc * 1 + 72)/2" | LC_ALL=C bc)

        fslmaths $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz -mul 0 -add 1 -roi $x_rhpc 1 $y_rhpc 1 $z_rhpc 1 0 1 ./overlays/${id}_righthpc -odt float
        fslmaths ./overlays/${id}_righthpc.nii.gz -kernel sphere $sphere_size -fmean ./overlays/${id}_righthpc_sphere -odt float
        fslmaths ./overlays/${id}_righthpc_sphere.nii.gz -bin ./overlays/${id}_righthpc_bin.nii.gz 
        fslmaths ./overlays/${id}_hpc_mask -add ./overlays/${id}_righthpc_bin.nii.gz ./overlays/${id}_hpc_mask
        rm ./overlays/${id}_righthpc.nii.gz ./overlays/${id}_righthpc_sphere.nii.gz ./overlays/${id}_righthpc_bin.nii.gz
    fi

    # combine across patients
    fslmaths ./overlays/amygdala_mask -add ./overlays/${id}_amygdala_mask ./overlays/amygdala_mask
    fslmaths ./overlays/amygdala_mask -bin ./overlays/amygdala_mask 
    fslmaths ./overlays/hpc_mask -add ./overlays/${id}_hpc_mask ./overlays/hpc_mask
    fslmaths ./overlays/hpc_mask -bin ./overlays/hpc_mask 
done < $position_file

# To reproduce the figure set fsleyes position voxel coordinates: 65 107 55

