function vgg16_layer_features = load_vgg16_layer_features(file_path, stim_used)
    vgg16_layer_features = readtable(file_path);
    % delete first column
    vgg16_layer_features(:,1) = [];
    % convert to array
    vgg16_layer_features = table2array(vgg16_layer_features);
    % exclude stimuli that are not used
    vgg16_layer_features(stim_used==0,:) = [];
    % replace NaNs with zeros
    vgg16_layer_features(isnan(vgg16_layer_features)) = 0;
end