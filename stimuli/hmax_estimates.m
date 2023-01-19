clc
clear
% path to HMAX functions
addpath('/Users/kamp/PhD/spikevar/repo/SpikeVar/toolboxes/hmaxMatlab');
% paths to stimuli
% Use minimally renamed verion since it makes sorting things easier
image_dir_base = "/Users/kamp/PhD/spikevar/stimuli_rn/";
save_dir_base = "/Users/kamp/PhD/spikevar/output/hmax_output/";
versions = ["newolddelay", "newolddelay2", "newolddelay3"];
total_img = 590;

%% set up categroy names etc
cat_dict = [{"houses", "landscapes", "mobility", "phones", "smallAnimal"},...
    {"fruit", "kids", "military", "space", "zzanimal"},...
    {"1cars", "2food","3people", "4spatial", "5animals"}];
cat_nums = 1:15;

%% Initialize S1 gabor filters and C1 parameters
fprintf('initializing S1 gabor filters\n');
orientations = [90 -45 0 45]; % 4 orientations for gabor filters
RFsizes      = 7:2:37;        % receptive field sizes
div          = 4:-.05:3.25;    % tuning parameters for the filters' "tightness"
[filterSizes,filters,c1OL,~] = initGabor(orientations,RFsizes,div);
fprintf('initializing C1 parameters\n')
c1Scale = 1:2:18; % defining 8 scale bands
c1Space = 8:2:22; % defining spatial pooling range for each scale band

%% Load the universal patch set.
fprintf('Loading the universal patch set\n')
load('universal_patch_set.mat','patches','patchSizes');
nPatchSizes = size(patchSizes,2);

%% loop across versions
HMAX_info = table(NaN(total_img,1), NaN(total_img,1), NaN(total_img,1), NaN(total_img,1), NaN(total_img,1), ... 
    'VariableNames', {'stimulus_code', 'version', 'category', 'c1_mean', 'c2_mean'});
HMAX_estim = struct('patch_wise_med_c1', NaN, 'patch_wise_med_c2', NaN, 'patch_wise_sd_c1', NaN, ...
    'patch_wise_sd_c2', NaN, 'ori_wise_med_c1', NaN, 'ori_wise_sd_c1', NaN, ...
    'ori_wise_med_c1_H', NaN, 'ori_wise_c1_H', NaN, 'img_H', NaN);

% initialize
patch_wise_med_c1 = []; patch_wise_med_c2 = [];
patch_wise_sd_c1 = []; patch_wise_sd_c2 = [];

n_patches = 8; n_orientations = 4;
ori_wise_med_c1 = NaN(total_img, n_patches, n_orientations); 
ori_wise_sd_c1 = NaN(total_img, n_patches, n_orientations);
ori_wise_med_c1_H = NaN(total_img, n_patches, n_orientations); 
ori_wise_c1_H = NaN(total_img, n_patches, n_orientations);

img_H = [];

fprintf('Starting HMAX\n')
idx = 1;
for ver = 1:length(versions)
    verdir = [image_dir_base + versions(ver) + '/'];
    % loop across categories
    for cat = 1:5
        catdir = [verdir + cat_dict(5*(ver-1)+cat) + '/'];
        % get list of images in the folder of interest
        pic_list = dir(catdir + '*.jpg');
        n_img = length(pic_list);
        
        % loop across pictures in folder, load and make grayscale
        for ipic = 1:n_img
            % get the unique stimulus number that is used later
            stim_code = 1000 + cat_nums(5*(ver-1)+cat)*100+ipic;
            % load image
            pic_path = [catdir + pic_list(ipic).name];
            picdat = load_pic(pic_path);
            
            % Extract layer responses
            [c2,c1] = extractC2forCell...
                (filters,filterSizes,c1Space,c1Scale,c1OL,patches,{picdat},nPatchSizes,patchSizes(1:3,:));
            
            % calculate responses for each patch
            patch_estim = get_patch_wise_estim(c1,c2);
            % store patch-wise C1 and C2 values
            patch_wise_med_c1 = [patch_wise_med_c1; patch_estim.c1_median];
            patch_wise_med_c2 = [patch_wise_med_c2; cell2mat(patch_estim.c2_median)];
            patch_wise_sd_c1 = [patch_wise_sd_c1; patch_estim.c1_sd];
            patch_wise_sd_c2 = [patch_wise_sd_c2; cell2mat(patch_estim.c2_sd)];
            % c1/c2 mean
            c1_mean = mean(patch_estim.c1_median);
            c2_mean = mean(cell2mat(patch_estim.c2_median));
            
            % get entropy and orientation-wise estimates
            ori_estim = get_ori_wise_estim(c1, c1Scale, orientations);
            % store orientation-wise estimates
            ori_wise_med_c1(idx,:,:) = ori_estim.c1_median; 
            ori_wise_sd_c1(idx,:,:) = ori_estim.c1_sd;
            ori_wise_med_c1_H(idx,:,:) = ori_estim.c1_median_entropy; 
            ori_wise_c1_H(idx,:,:) =  ori_estim.c1_entropy;
            
            % entropy of image 
            entropy_img = entropy(picdat);
           
            % store info to table
            HMAX_info(idx,:) = {stim_code, ver, cat, ... 
                c1_mean, c2_mean};
            
            idx = idx + 1;
            clear c1
            clear c2 
            
        end
    end
    
end

% store all data in one table
writetable(HMAX_info, [save_dir_base + 'HMAX_info.csv'])

save([save_dir_base + 'HMAX_estim.mat'], 'HMAX_info',...
    'patch_wise_med_c1', 'patch_wise_med_c2', 'patch_wise_sd_c1', 'patch_wise_sd_c2',...
    'ori_wise_med_c1', 'ori_wise_sd_c1', 'ori_wise_med_c1_H', 'ori_wise_c1_H',...
    'img_H');