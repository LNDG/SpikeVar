"""
Created on Tue Apr 12 09:40:41 2022
@author: waschke & kamp
Here we calculate visual complexity by calculating the mean value of layer 4 output in VGG16
"""

import tensorflow as tf
import numpy as np
from tensorflow import keras
import os
import pandas as pd
import utils

# set up folder structure
base_dir = "/Users/kamp/PhD/spikevar"
image_dir_base = os.path.join(base_dir, "stimuli_rn")
output_dir = os.path.join(base_dir, "output", "vgg16_output")

# versions and subcategories
version_dict = {
    "newolddelay": {"houses": 1, "landscapes": 2, "mobility": 3, "phones": 4, "smallAnimal": 5},
    "newolddelay2": {"fruit": 6, "kids": 7, "military": 8, "space": 9, "zzanimal": 10},
    "newolddelay3": {"1cars":11, "2food": 12,"3people":13, "4spatial": 14, "5animals": 15}
}

# get image paths
img_dict, img_ls = utils.get_image_paths(image_dir_base, version_dict)
num_images = len(img_ls)
stimulus_codes_df = []


# build model
model = tf.keras.applications.VGG16(
    include_top=True,
    weights="imagenet",
    input_tensor=None,
    input_shape=None,
    pooling=None,
#   classes=1000,
    classifier_activation="softmax",
)

# instantiating a activation model from a input tensor and all imtermidate output tensor
layer_outputs = [layer.output for layer in model.layers[:]]
activation_model = keras.models.Model(inputs=model.input, outputs=layer_outputs)

# layers of interest are the max-pooling layers
layer_names = ['first', 'second', 'third', 'fourth', 'fifth', 'final']
layer_idx = {'first': 3, 'second': 6, 'third': 10, 'fourth':14, 'fifth': 18, 'final': 21}
lw_feature_count = [64, 128, 256, 512, 512, 4096]
layer_sizes = [112, 56, 28, 14, 7, 1]

# set up arrays for feature maps 
layer_activations = {layer_names[i]: np.empty([0, layer_sizes[i], layer_sizes[i], lw_feature_count[i]]) 
                for i in range(5)}
layer_activations.update(final = np.empty([0, 4096]))
layer_activations.update(all_vers = np.empty([0]))

# set up arrays for complexity estimates
complexity_dict = {}
complexity_dict.update({f'{layer}_layer_mean': np.empty([num_images,]) for layer in layer_names})
complexity_dict.update({f'{layer}_layer_sd': np.empty([num_images,]) for layer in layer_names})

# set up arrays for non-zero complexity
complexity_dict.update({f'{layer}_layer_non-zero_mean': np.empty([num_images,]) for layer in layer_names})
complexity_dict.update({f'{layer}_layer_non-zero_sd': np.empty([num_images,]) for layer in layer_names})

# set up array for all features
lw_feature_complexity_mean = {layer_names[i]: np.empty([num_images, lw_feature_count[i]]) for i in range(5)}
lw_feature_complexity_sd = {layer_names[i]: np.empty([num_images, lw_feature_count[i]]) for i in range(5)}
lw_feature_complexity_sum = {layer_names[i]: np.empty([num_images, lw_feature_count[i]]) for i in range(5)}

# non-zero
lw_feature_nz_complexity_mean = {layer_names[i]: np.empty([num_images, lw_feature_count[i]]) for i in range(5)}
lw_feature_nz_complexity_sd = {layer_names[i]: np.empty([num_images, lw_feature_count[i]]) for i in range(5)}

# get the layerwise feature maps for each iamge
for version in version_dict.keys():
    for category in version_dict[version]:
        img_ls = img_dict[version][category]
        batch_size = len(img_ls)
        imgs_batch_tensor = np.empty((batch_size, 224, 224, 3))  
        category_code = version_dict[version][category]      
        for i in range(batch_size):
            # get stimulus code
            image_number = i + 1
            image_path = img_ls[i]
            stimulus_code = utils.get_stimulus_code(image_number, category_code, category, version, image_path)
            stimulus_codes_df.append(stimulus_code)
            # load image and check size compatibility
            img_tensor_rgb = utils.load_image_tensor(image_path)
            assert img_tensor_rgb.size == (224, 224), f"img_tensor.size:{img_tensor_rgb.size}"
            imgs_batch_tensor[i, :, :, :] = img_tensor_rgb
        
        # predict activation maps for all images in this batch
        activations = activation_model.predict(imgs_batch_tensor)
        
        # key part: get and store feature maps of layers of interest 
        for layer in layer_names:
            layer_activations[layer] = np.append(layer_activations[layer], activations[layer_idx[layer]], axis=0)     
        print(f"Done with estimates for version: {version} and category: {category}")

stimulus_codes_df = pd.DataFrame(stimulus_codes_df)

for i in range(num_images):
    for j, layer in enumerate(layer_names):
        img_lw_activations = layer_activations[layer][i,...]
        # get mean/sd of each layer
        complexity_dict[f'{layer}_layer_mean'][i] = np.mean(img_lw_activations)        
        complexity_dict[f'{layer}_layer_sd'][i] = np.std(img_lw_activations)
        # due to sparsity, get the mean only across non-zerop entries
        complexity_dict[f'{layer}_layer_non-zero_mean'][i] = np.mean(img_lw_activations[np.nonzero(img_lw_activations)])
        complexity_dict[f'{layer}_layer_non-zero_sd'][i] = np.std(img_lw_activations[np.nonzero(img_lw_activations)])
        
        if layer == 'final': 
            continue

        for n_feature in range(lw_feature_count[j]):
            feature_activations = img_lw_activations[..., n_feature]
            # get mean/sd/sum of each feature in each layer
            lw_feature_complexity_mean[layer][i, n_feature] = np.nanmean(feature_activations)
            lw_feature_complexity_sd[layer][i, n_feature] = np.std(feature_activations)
            lw_feature_complexity_sum[layer][i, n_feature] = np.sum(feature_activations)
            # same for non-zero entries
            lw_feature_nz_complexity_mean[layer][i, n_feature] = np.nanmean(feature_activations[np.nonzero(feature_activations)])
            lw_feature_nz_complexity_sd[layer][i, n_feature] = np.std(feature_activations[np.nonzero(feature_activations)])
             
# save layerwise mean/sd/sum to csv
stim_dat = pd.DataFrame(complexity_dict)
stim_dat = stimulus_codes_df.join(stim_dat)
stim_dat.to_csv(output_dir + "SpikeVar_VGG16_Complexity.csv")
                              
# save layer-wise feature maps
for layer_name, feature_count in zip(layer_names[:-1], lw_feature_count[:-1]): 
    col_names = [f"Feature_{i:04d}" for i in range(feature_count)]
    pd.DataFrame(lw_feature_complexity_mean[layer_name], columns=col_names).to_csv(os.path.join(output_dir, f"SpikeVar_VGG16_{layer_name}_layer_features_mean.csv"))
    pd.DataFrame(lw_feature_complexity_sd[layer_name], columns=col_names).to_csv(os.path.join(output_dir, f"SpikeVar_VGG16_{layer_name}_layer_features_sd.csv"))
    pd.DataFrame(lw_feature_complexity_sum[layer_name], columns=col_names).to_csv(os.path.join(output_dir, f"SpikeVar_VGG16_{layer_name}_layer_features_sum.csv"))
    # non-zero estimates
    pd.DataFrame(lw_feature_nz_complexity_mean[layer_name], columns=col_names).to_csv(os.path.join(output_dir, f"SpikeVar_VGG16_{layer_name}_layer_non-zero_features_mean.csv"))
    pd.DataFrame(lw_feature_nz_complexity_sd[layer_name], columns=col_names).to_csv(os.path.join(output_dir, f"SpikeVar_VGG16_{layer_name}_layer_non-zero_features_sd.csv"))
    # save the complete feature space per layer
    np.save(os.path.join(output_dir, f"SpikeVar_VGG16_{layer_name}_layer_all_features"), layer_activations[layer_name])
