import os 
import pandas as pd

def get_image_paths(root, versions, cat_names):
    """
    param root: base directory of the stimuli files
    param versions: list of all the versions
    param cat_names: list of categories separated for each version
    """ 
    image_dict = {versions[i]: {cat_names[i][j]: sorted([os.path.join(root, versions[i], cat_names[i][j], image) 
        for image in os.listdir(os.path.join(root, versions[i], cat_names[i][j])) if image.endswith('.jpg')]) 
        for j in range(5)} 
        for i in range(3)}
    
    return image_dict

def get_stimulus_codes(image_dict):
    cat_dict = {
        "newolddelay": {"houses": 1, "landscapes": 2, "mobility": 3, "phones": 4, "smallAnimal": 5},
        "newolddelay2": {"fruit":6, "kids": 7, "military": 8, "space": 9, "zzanimal": 10},
        "newolddelay3": {"1cars":11, "2food": 12,"3people":13, "4spatial": 14, "5animals": 15}
        }

    code_df = []
    for version in image_dict.keys():
        for category in image_dict[version]:
            cat_code = cat_dict[version][category]
            for i, image in enumerate(image_dict[version][category]):
                image_number = i + 1 
                stim_code =  1000 + (cat_code * 100) + image_number
                code_df.append({'stimulus_code': stim_code, 
                    'category': category, 
                    'version': version, 
                    'image_number': image_number, 
                    'image_path': image})
    return pd.DataFrame(code_df)

def load_image_tensor(image_path):
    img_tensor = Image.open(image_path)  
    img_tensor_resize = img_tensor.resize((224,224))
    img_tensor_rgb = img_tensor_resize.convert('RGB')
    return img_tensor_rgb