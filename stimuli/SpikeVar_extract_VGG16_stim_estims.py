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
import utils

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

img_dict, img_ls = utils.get_image_paths(image_dir_base, versions, cat_names)
num_images = len(img_ls)
stimulus_codes_df = []
stimulus_codes_ls = []

layer_names = ['first', 'second', 'third', 'fourth', 'fifth', 'final']
layer_idx = {'first': 3, 'second': 6, 'third': 10, 'fourth':14, 'fifth': 18, 'final': 21}
lw_feature_count = [64, 128, 256, 512, 512, 4096]
layer_sizes = [112, 56, 28, 14, 7, 1]

# set up arrays for feature maps 
layer_activations = {layer_names[i]: np.empty([0, layer_sizes[i], layer_sizes[i], lw_feature_count[i]] 
                for i in range(5))}
layer_activations.update(final = np.empty([0, 4096]))
layer_activations.update(all_vers = np.empty([0]))

# set up arrays for complexity estimates
complexity_dict = {}
complexity_dict.update({f'{layer}_layer_mean': np.empty([590,1]) for layer in layer_names})
complexity_dict.update({f'{layer}_layer_sd': np.empty([590,1]) for layer in layer_names})

# set up arrays for non-zero complexity
complexity_dict.update({f'{layer}_layer_non-zero_mean': np.empty([590,1]) for layer in layer_names})
complexity_dict.update({f'{layer}_layer_non-zero_sd': np.empty([590,1]) for layer in layer_names})

# set up array for all features
lw_feature_complexity_mean = {layer_names[i]: np.empty([590, lw_feature_count[i]]) for i in range(5)}
lw_feature_complexity_sd = {layer_names[i]: np.empty([590, lw_feature_count[i]]) for i in range(5)}
lw_feature_complexity_sum = {layer_names[i]: np.empty([590, lw_feature_count[i]]) for i in range(5)}

# non-zero
lw_feature_nz_complexity_mean = {layer_names[i]: np.empty([590, lw_feature_count[i]]) for i in range(5)}
lw_feature_nz_complexity_sd = {layer_names[i]: np.empty([590, lw_feature_count[i]]) for i in range(5)}

# get the layerwise feature maps for each iamge
for ver in versions:
    for cat in range(len(cat_names[ver])):
        img_ls = img_dict[ver][cat]
        batch_size = len(img_ls)
        imgs_batch_tensor = np.empty((batch_size, 224, 224, 3))  
        cat_code = cat_dict[version][category]      
        for i in range(batch_size):
            # get stimulus code
            image_number = i + 1 
            stim_code =  1000 + (cat_code * 100) + image_number
            image_path = img_ls[i]
            stimulus_codes_df.append({'stimulus_code': stim_code, 
                'category': category, 
                'version': version, 
                'image_number': image_number, 
                'image_path': image})
            img_tensor_rgb = utils.load_image_tensor(image_path)
            # check size compatibility
            assert img_tensor_rgb.shape == (224, 224, 3), f"img_tensor.size:{img_tensor_rgb.size}"
            imgs_batch_tensor[i, :, :, :] = img_tensor_rgb
        
        # predict activation maps for all images in this batch
        activations = activation_model.predict(imgs_batch_tensor)
        
        # key part: get and store feature maps of layers of interest 
        for layer in layer_names:
            np.append(layer_activations[layer], (activations[layer_idx[layer]]), axis=0)
              
        print(f"Done with estimates for version {ver}, category {cat + 1}")

for i in range(num_images):
    for j, layer in enumerate(layer_names):
        img_lw_activations = layer_activations[layer][i,...]
        # get mean/sd of each layer
        complexity_dict[f'{layer}_layer_mean'][i] = np.mean(img_lw_activations)         
        complexity_dict[f'{layer}_layer_sd'][i] = np.std(img_lw_activations)
        # due to sparsity, get the mean only across non-zerop entries
        complexity_dict[f'{layer}_layer_non-zero_mean'][i] = np.mean(img_lw_activations[np.nonzero(img_lw_activations)])
        complexity_dict[f'{layer}_layer_non-zero_sd'][i] = np.std(img_lw_activations[np.nonzero(img_lw_activations)])
        
        for n_feature in range(lw_feature_count[j]):
            feature_activations = img_lw_activations[i,:,:, n_feature]
            # get mean/sd/sum of each feature in each layer
            lw_feature_complexity_mean[layer][i, n_feature] = np.nanmean(feature_activations)
            lw_feature_complexity_sd[layer][i, n_feature] = np.std(feature_activations)
            lw_feature_complexity_sum[layer][i, n_feature] = np.sum(feature_activations)
            # same for non-zero entries
            lw_feature_nz_complexity_mean[layer][i, n_feature] = np.nanmean(feature_activations[np.nonzero(feature_activations)])
            lw_feature_nz_complexity_sd[layer][i, n_feature] = np.std(feature_activations[np.nonzero(feature_activations)])
             
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
