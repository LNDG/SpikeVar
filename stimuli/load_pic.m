function picdat = load_pic(path_to_pic)
% load image, check if grayscale already
    temp_img = imread(path_to_pic);
    sz = size(temp_img);
    % either use gray scale image or transform RGB to grayscale
    % in either case, resize to 224*224
    % this is done to keep things comparable across images and
    % across methods since VGG16 only eats 224*224 images rn
    if size(sz)<3
        picdat = double(imresize(temp_img,[224, 224]));
    else
        picdat = double(imresize(rgb2gray(temp_img),[224,224]));
    end
end
            