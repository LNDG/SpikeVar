import os 
import pandas as pd
from PIL import Image

def get_image_paths(root, versions_dict):
    """
    param root: base directory of the stimuli files
    param versions: list of all the versions
    param cat_names: list of categories separated for each version
    """ 
    img_dict = {version: {} for version in versions_dict.keys()}
    for version in versions_dict.keys():
        for category in versions_dict[version].keys():
            img_dict[version].update({category: sorted([os.path.join(root, version, category, image) 
                for image in os.listdir(os.path.join(root, version, category)) if image.endswith('.jpg')])})
    
    img_list = [image_path for version in img_dict.keys() 
                     for category in img_dict[version].keys() 
                     for image_path in img_dict[version][category]]
    
    return img_dict, img_list

def get_stimulus_code(image_number, category_code, category, version, image_path):
    version_numbers = {"newolddelay":1, "newolddelay2": 2, "newolddelay3": 3}
    stim_code =  1000 + (category_code * 100) + image_number
    stimulus_dict = {'stimulus_code': stim_code, 
        'category': category, 
        'version': version, 
        'version_number': version_numbers[version],
        'image_number': image_number, 
        'image_path': image_path}
    return stimulus_dict

def load_image_tensor(image_path):
    img_tensor = Image.open(image_path)  
    img_tensor_resize = img_tensor.resize((224,224))
    img_tensor_rgb = img_tensor_resize.convert('RGB')
    return img_tensor_rgb