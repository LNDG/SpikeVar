#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 09:40:41 2022

@author: waschke & kamp

Here we Calculate visual complexity by calculating the mean value of layer 4 output in VGG16
key part in line 54-55
"""

import tensorflow as tf
import numpy as np
from PIL import Image
from tensorflow import keras
import matplotlib.pyplot as plt
import csv
import cv2
import glob
import os
import pandas as pd

# build model
model = tf.keras.applications.VGG16(
    include_top=True,
    weights="imagenet",
    input_tensor=None,
    input_shape=None,
    pooling=None,
#    classes=1000,
    classifier_activation="softmax",
)
model.summary()

# instantiating a activation model from a input tensor and all imtermidate output tensor
layer_outputs = [layer.output for layer in model.layers[:]]
activation_model = keras.models.Model(inputs=model.input, outputs=layer_outputs)

#set up folder structure and remember:
    # we have three differtn task versions (newolddelay, newolddelay2, newolddelay3)
    # every version has its own subfolders for stimulus categories


# Visual complexity calculation
image_dir_base = "/Users/waschke/OneDrive/Old data and code/Old data and code/RecogMemory_MTL_release_v2/Stimuli_rn/"
save_dir_base = "/Users/waschke/OneDrive/Old data and code/VGG16_output/"
versions = ["newolddelay", "newolddelay2", "newolddelay3"]

cat_dict = [
    {"houses": 1, "landscapes": 2, "mobility": 3, "phones": 4, "smallAnimal": 5},
    {"fruit":6, "kids": 7, "military": 8, "space": 9, "zzanimal": 10},
    {"1cars":11, "2food": 12,"3people":13, "4spatial": 14, "5animals": 15}
        ]


cat_names = [
         ["houses", "landscapes", "mobility", "phones", "smallAnimal"],
         ["fruit", "kids", "military", "space", "zzanimal"],
        ["1cars", "2food","3people", "4spatial", "5animals"]
        ]

cat_list = []
stim_ID_list = []
version_list = []
img_list = []

layers = ['first', 'second', 'third', 'fourth', 'fifth', 'final']
lw_complexity_mean = {layer: np.empty([590,1]) for layer in layers}
lw_complexity_sd = {layer: np.empty([590,1]) for layer in layers}

# non-zero
lw_nz_complexity_mean = {layer: np.empty([590,1]) for layer in layers} # final complexity_list
lw_nz_complexity_sd = {layer: np.empty([590,1]) for layer in layers}

lw_feature_count = [64, 128, 256, 512, 512, 4096]
layer_sizes = [112, 56, 28, 14, 7, 1]

layer_activations = {layer[i]: np.empty([0, layer_sizes[i], layer_sizes[i], lw_feature_count[i]] 
                for i in range(5))}
layer_activations.update(final = np.empty([0, 4096]))
layer_activations.update(all_vers=np.empty([0]))

lw_complexity = {layer[i]: np.empty([590, lw_feature_count[i]]) for i in range(5)}


first_lw_wz_complexity_array = np.empty([590,64])
second_lw_wz_complexity_array = np.empty([590,128])
third_lw_wz_complexity_array = np.empty([590,256])
fourth_lw_wz_complexity_array = np.empty([590,512])
fifth_lw_wz_complexity_array = np.empty([590,512])

first_lw_wz_sd_array = np.empty([590,64])
second_lw_wz_sd_array = np.empty([590,128])
third_lw_wz_sd_array = np.empty([590,256])
fourth_lw_wz_sd_array = np.empty([590,512])
fifth_lw_wz_sd_array = np.empty([590,512])

first_lw_sd_array = np.empty([590,64])
second_lw_sd_array = np.empty([590,128])
third_lw_sd_array = np.empty([590,256])
fourth_lw_sd_array = np.empty([590,512])
fifth_lw_sd_array = np.empty([590,512])

first_lw_sum_array = np.empty([590,64])
second_lw_sum_array = np.empty([590,128])
third_lw_sum_array = np.empty([590,256])
fourth_lw_sum_array = np.empty([590,512])
fifth_lw_sum_array = np.empty([590,512])


# get number of files for each task variant and picture category
for ver in range(len(versions)):
    dd = image_dir_base + versions[ver]
    subcat_list = [ f.path for f in os.scandir(dd) if f.is_dir() ]
    cur_catnames = cat_names[ver]
    for scat in range(len(cat_names[ver])):
        #img_list = os.listdir(subcat_list[scat])    
        # get index of image category which is not identical with the running index
        #temp_p_dis = os.path.split(subcat_list[scat])
        catname = cur_catnames[scat]
        
        
        img_list_path = sorted(glob.glob(dd+"/"+catname + "/*.jpg", recursive = False))
        print("Working on Version {}".format(ver+1))
        img_ls = []
        for img in range(len(img_list_path)):
            img_ls.append(os.path.basename(img_list_path[img]))
        batch_size = len(img_ls)
        imgs_batch_tensor = np.empty((batch_size, 224, 224, 3))        
        for i in range(batch_size):
            img_tensor = Image.open(img_list_path[i])  
            img_tensor_resize = img_tensor.resize((224,224))
            img_tensor_rgb = img_tensor_resize.convert('RGB')
            # check size compatibility
            try:
                imgs_batch_tensor[i, :, :, :] = img_tensor_rgb
            except:
                print(f"img_tensor.size:{img_tensor.size}")
                raise NotImplementedError
                
            # determine and set up stimulus code
            cat_code = cat_dict[ver][catname]
            stim_code =  1000+(cat_code*100)+i+1
            stim_ID_list.append(stim_code)
            cat_list.append(catname)    
            version_list.append(ver+1)
            img_list.append(i+1)
            all_vers = np.append(all_vers, np.array(ver+1))
            
        activations = activation_model.predict(imgs_batch_tensor)
        final_layer_activation = activations[21]
        fifth_layer_activation = activations[18]
        fourth_layer_activation = activations[14]  # key part: get layer 4 feature maps, and its mean value can represent complexity well
        third_layer_activation = activations[10]
        second_layer_activation = activations[6]
        first_layer_activation = activations[3]
        
        # append stuff
        all_first_layer_activations  =np.append(all_first_layer_activations, (first_layer_activation), axis = 0)
        all_second_layer_activations  =np.append(all_second_layer_activations, (second_layer_activation), axis = 0)
        all_third_layer_activations  =np.append(all_third_layer_activations, (third_layer_activation), axis = 0)
        all_fourth_layer_activations  =np.append(all_fourth_layer_activations, (fourth_layer_activation), axis = 0)
        all_fifth_layer_activations  =np.append(all_fifth_layer_activations, (fifth_layer_activation), axis = 0)
        all_final_layer_activations  =np.append(all_final_layer_activations, (final_layer_activation), axis = 0)
        
      
        print("Done with estimates for version {}, category {}".format(ver+1,scat+1))
        
for i in range(all_second_layer_activations.shape[0]):
    final_complexity_list[i] = ([np.mean(all_final_layer_activations[i, :]) ]) 
    fifth_complexity_list[i] =([np.mean(all_fifth_layer_activations[i, :, :, :]) ]) 
    fourth_complexity_list[i] = ([np.mean(all_fourth_layer_activations[i, :, :, :]) ]) 
    third_complexity_list[i]=([np.mean(all_third_layer_activations[i, :, :, :]) ] )
    second_complexity_list[i]=([np.mean(all_second_layer_activations[i, :, :, :])]) 
    first_complexity_list[i]=([np.mean(all_first_layer_activations[i, :, :, :])]) 
        
    # due to sparsity, why not get the mean only across non-zerop entries?
    final_nz_complexity_list[i] = ( [np.mean(np.nonzero(all_final_layer_activations[i, :].flatten()))])
    fifth_nz_complexity_list[i] = ([np.mean(np.nonzero(all_fifth_layer_activations[i, :, :, :].flatten()))]) 
    fourth_nz_complexity_list[i] = ([np.mean(np.nonzero(all_fourth_layer_activations[i, :, :, :].flatten()))]) 
    third_nz_complexity_list[i] = ( [np.mean(np.nonzero(all_third_layer_activations[i, :, :, :].flatten()))] )
    second_nz_complexity_list[i] = ([np.mean(np.nonzero(all_second_layer_activations[i, :, :, :].flatten()))] )
    first_nz_complexity_list[i] = ([np.mean(np.nonzero(all_first_layer_activations[i, :, :, :].flatten()))] )
    # can we get the sd on each layer?
    final_complexity_sd_list[i] = ([np.std(all_final_layer_activations[i, :])] )
    fifth_complexity_sd_list[i] = ( [np.std(all_fifth_layer_activations[i, :, :, :])] )
    fourth_complexity_sd_list[i] = ([np.std(all_fourth_layer_activations[i, :, :, :])] )
    third_complexity_sd_list[i] = ( [np.std(all_third_layer_activations[i, :, :, :])] )
    second_complexity_sd_list[i] = ([np.std(all_second_layer_activations[i, :, :, :])])
    first_complexity_sd_list[i] = ([np.std(all_first_layer_activations[i, :, :, :])])
        
            # sd for nonzero elements?
    final_nzsd_complexity_list[i] = ( [np.std(np.nonzero(all_final_layer_activations[i, :].flatten()))] )
    fifth_nzsd_complexity_list[i] = ( [np.std(np.nonzero(all_fifth_layer_activations[i, :, :, :].flatten())) ] )
    fourth_nzsd_complexity_list[i] = ( [np.std(np.nonzero(all_fourth_layer_activations[i, :, :, :].flatten())) ] )
    third_nzsd_complexity_list[i] = ([np.std(np.nonzero(all_third_layer_activations[i, :, :, :].flatten())) ] )
    second_nzsd_complexity_list[i] = ([np.std(np.nonzero(all_second_layer_activations[i, :, :, :].flatten())) ] )
    first_nzsd_complexity_list[i] = ([np.std(np.nonzero(all_first_layer_activations[i, :, :, :].flatten())) ] )
        
            # one more thing to try: reduce feature space by running a PCA on the layer-specific output
            # problem: many zero entries and hence hard to normalize

            # before that: average (non-zero) layer outputs within features but across image dimensions
            
            
    for nf1 in range(first_layer_activation.shape[-1]):
      first_lw_complexity_array[i, nf1]  = np.nanmean(np.nonzero(all_first_layer_activations[i, :, :, nf1].flatten())) 
      first_lw_wz_complexity_array[i, nf1]  = np.nanmean((all_first_layer_activations[i, :, :, nf1].flatten())) 
      # also get the sd for each feature as a proxy of sparsity
      first_lw_wz_sd_array[i, nf1]  = np.std((all_first_layer_activations[i, :, :, nf1].flatten())) 
      first_lw_sd_array[i, nf1]  = np.std(np.nonzero(all_first_layer_activations[i, :, :, nf1].flatten())) 
      # one more thing: get the usm of non-zero elements
      first_lw_sum_array[i, nf1]  = np.sum(np.nonzero(all_first_layer_activations[i, :, :, nf1].flatten())) 
    for nf2 in range(second_layer_activation.shape[-1]):
        second_lw_complexity_array[i, nf2]  = np.nanmean(np.nonzero(all_second_layer_activations[i, :, :, nf2].flatten())) 
        second_lw_wz_complexity_array[i, nf2]  = np.nanmean((all_second_layer_activations[i, :, :, nf2].flatten())) 
        # also get the sd for each feature as a proxy of sparsity
        second_lw_wz_sd_array[i, nf2]  = np.std((all_second_layer_activations[i, :, :, nf2].flatten())) 
        second_lw_sd_array[i, nf2]  = np.std(np.nonzero(all_second_layer_activations[i, :, :, nf2].flatten())) 
        # one more thing: get the usm of non-zero elements
        second_lw_sum_array[i, nf2]  = np.sum(np.nonzero(all_second_layer_activations[i, :, :, nf2].flatten())) 
        
    for nf3 in range(third_layer_activation.shape[-1]):
        third_lw_complexity_array[i, nf3]  = np.nanmean(np.nonzero(all_third_layer_activations[i, :, :, nf3].flatten()))  
        third_lw_wz_complexity_array[i, nf3]  = np.nanmean((all_third_layer_activations[i, :, :, nf3].flatten()))  
        third_lw_wz_sd_array[i, nf3]  = np.std((all_third_layer_activations[i, :, :, nf3].flatten())) 
        third_lw_sd_array[i, nf3]  = np.std(np.nonzero(all_third_layer_activations[i, :, :, nf3].flatten()))  
        third_lw_sum_array[i, nf3]  = np.sum(np.nonzero(all_third_layer_activations[i, :, :, nf3].flatten()))  
    for nf4 in range(fourth_layer_activation.shape[-1]):
        fourth_lw_complexity_array[i, nf4]  = np.nanmean(np.nonzero(all_fourth_layer_activations[i, :, :, nf4].flatten()))  
        fourth_lw_wz_complexity_array[i, nf4]  = np.nanmean((all_fourth_layer_activations[i, :, :, nf4].flatten()))  
        fourth_lw_wz_sd_array[i, nf4]  = np.std((all_fourth_layer_activations[i, :, :, nf4].flatten()))  
        fourth_lw_sd_array[i, nf4]  = np.std(np.nonzero(all_fourth_layer_activations[i, :, :, nf4].flatten()))  
        fourth_lw_sum_array[i, nf4]  = np.sum(np.nonzero(all_fourth_layer_activations[i, :, :, nf4].flatten()))  
    for nf5 in range(fifth_layer_activation.shape[-1]):
        fifth_lw_complexity_array[i, nf5]  = np.nanmean(np.nonzero(all_fifth_layer_activations[i, :, :, nf5].flatten()))  
        fifth_lw_wz_complexity_array[i, nf5]  = np.nanmean((all_fifth_layer_activations[i, :, :, nf5].flatten()))  
        fifth_lw_wz_sd_array[i, nf5]  = np.std((all_fifth_layer_activations[i, :, :, nf5].flatten()))  
        fifth_lw_sd_array[i, nf5]  = np.std(np.nonzero(all_fifth_layer_activations[i, :, :, nf5].flatten()))    
        fifth_lw_sum_array[i, nf5]  = np.sum(np.nonzero(all_fifth_layer_activations[i, :, :, nf5].flatten()))  

       

stim_Dat = np.column_stack((second_complexity_list, third_complexity_list,fourth_complexity_list, fifth_complexity_list,
                            final_complexity_list, 
                            second_nz_complexity_list, third_nz_complexity_list,fourth_nz_complexity_list, fifth_nz_complexity_list,
                            final_nz_complexity_list,
                            second_complexity_sd_list, third_complexity_sd_list,fourth_complexity_sd_list, 
                            fifth_complexity_sd_list, final_complexity_sd_list, 
                            second_nzsd_complexity_list, third_nzsd_complexity_list,fourth_nzsd_complexity_list, fifth_nzsd_complexity_list,
                            final_nzsd_complexity_list,
                            stim_ID_list, cat_list,version_list,img_list))

# write results to csv
pd.DataFrame(stim_Dat).to_csv(save_dir_base + "SpikeVar_VGG16_Complexity.csv", header = ['Complexity_2', 'Complexity_3', 'Complexity_4',
                                                                                         'Complexity_5', 'Complexity_final',
                                                                                         'NZ_Complexity_2', 'NZ_Complexity_3', 'NZ_Complexity_4',
                                                                                         'NZ_Complexity_5', 'NZ_Complexity_final',
                                                                                         'Complexity_sd_2', 'Complexity_sd_3',
                                                                                         'Complexity_sd_4', 'Complexity_sd_5',
                                                                                         'Complexity_sd_final',
                                                                                         'Complexity_nzsd_2', 'Complexity_nzsd_3',
                                                                                         'Complexity_nzsd_4', 'Complexity_nzsd_5',
                                                                                         'Complexity_nzsd_final',
                                                                                         'StimulusCode', 'StimulusCategory',
                                                                                         'Version', 'Stimulus'])
                              
# generate names so we can save things using pandas

num_cols_1st = first_lw_complexity_array.shape[1]
rng = range(num_cols_1st)
new_cols = ['Feature_' + str(i+1) for i in range(num_cols_1st)]
pd.DataFrame(first_lw_complexity_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay1_features.csv", header = new_cols)
pd.DataFrame(first_lw_wz_complexity_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay1_wz_features.csv", header = new_cols)
pd.DataFrame(first_lw_wz_sd_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay1_wz_feature_sd.csv", header = new_cols)
pd.DataFrame(first_lw_sd_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay1_feature_sd.csv", header = new_cols)
pd.DataFrame(first_lw_sum_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay1_feature_sum.csv", header = new_cols)



num_cols_2nd = second_lw_complexity_array.shape[1]
rng = range(num_cols_2nd)
new_cols = ['Feature_' + str(i+1) for i in range(num_cols_2nd)]
pd.DataFrame(second_lw_complexity_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay2_features.csv", header = new_cols)
pd.DataFrame(second_lw_wz_complexity_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay2_wz_features.csv", header = new_cols)
pd.DataFrame(second_lw_wz_sd_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay2_wz_feature_sd.csv", header = new_cols)
pd.DataFrame(second_lw_sd_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay2_feature_sd.csv", header = new_cols)
pd.DataFrame(second_lw_sum_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay2_feature_sum.csv", header = new_cols)

num_cols_3rd = third_lw_complexity_array.shape[1]
rng = range(num_cols_3rd)
new_cols = ['Feature_' + str(i+1) for i in range(num_cols_3rd)]
pd.DataFrame(third_lw_complexity_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay3_features.csv", header = new_cols)
pd.DataFrame(third_lw_wz_complexity_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay3_wz_features.csv", header = new_cols)
pd.DataFrame(third_lw_wz_sd_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay3_wz_feature_sd.csv", header = new_cols)
pd.DataFrame(third_lw_sd_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay3_feature_sd.csv", header = new_cols)
pd.DataFrame(third_lw_sum_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay3_feature_sum.csv", header = new_cols)

num_cols_4th = fourth_lw_complexity_array.shape[1]
rng = range(num_cols_4th)
new_cols = ['Feature_' + str(i+1) for i in range(num_cols_4th)]
pd.DataFrame(fourth_lw_complexity_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay4_features.csv", header = new_cols)
pd.DataFrame(fourth_lw_wz_complexity_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay4_wz_features.csv", header = new_cols)
pd.DataFrame(fourth_lw_wz_sd_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay4_wz_feature_sd.csv", header = new_cols)
pd.DataFrame(fourth_lw_sd_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay4_feature_sd.csv", header = new_cols)
pd.DataFrame(fourth_lw_sum_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay4_feature_sum.csv", header = new_cols)

num_cols_5th = fifth_lw_complexity_array.shape[1]
rng = range(num_cols_5th)
new_cols = ['Feature_' + str(i+1) for i in range(num_cols_5th)]
pd.DataFrame(fifth_lw_complexity_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay5_features.csv", header = new_cols)
pd.DataFrame(fifth_lw_wz_complexity_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay5_wz_features.csv", header = new_cols)
pd.DataFrame(fifth_lw_wz_sd_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay5_wz_feature_sd.csv", header = new_cols)
pd.DataFrame(fifth_lw_sd_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay5_feature_sd.csv", header = new_cols)
pd.DataFrame(fifth_lw_sum_array).to_csv(save_dir_base + "SpikeVar_VGG16_lay5_feature_sum.csv", header = new_cols)


num_cols_final = all_final_layer_activations.shape[1]
rng = range(num_cols_final)
new_cols = ['Feature_' + str(i+1) for i in range(num_cols_final)]
pd.DataFrame(all_final_layer_activations).to_csv(save_dir_base + "SpikeVar_VGG16_lay_final_features.csv", header = new_cols)
# save task variant info as well
pd.DataFrame(all_vers).to_csv(save_dir_base + "SpikeVar_stim_task_variants.csv", header = ['Task_Variant'])


# save the whole feature space per layer
np.save(save_dir_base + "SpikeVar_VGG16_lay1_all_features", all_first_layer_activations)
np.save(save_dir_base + "SpikeVar_VGG16_lay2_all_features", all_second_layer_activations)
np.save(save_dir_base + "SpikeVar_VGG16_lay3_all_features", all_third_layer_activations)
np.save(save_dir_base + "SpikeVar_VGG16_lay4_all_features", all_fourth_layer_activations)
np.save(save_dir_base + "SpikeVar_VGG16_lay5_all_features", all_fifth_layer_activations)
