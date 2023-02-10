clc
clear
% path to HMAX functions
addpath('/Users/kamp/PhD/spikevar/repo/SpikeVar/toolboxes/hmaxMatlab');
% paths to stimuli
% Use minimally renamed verion since it makes sorting things easier
image_dir_base = "/Users/kamp/PhD/spikevar/stimuli_rn/";
save_dir_base = "/Users/kamp/PhD/spikevar/repo/SpikeVar/output/hmax/";
versions = ["newolddelay", "newolddelay2", "newolddelay3"];

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
stim_codes = []; pic_paths = [];
stim_versions = []; stim_categories = [];
c1_mean = []; c2_mean = [];
patch_wise_med_c1 = []; patch_wise_med_c2 = [];
patch_wise_sd_c1 = []; patch_wise_sd_c2 = [];
patch_wise_sum_c1 = []; patch_wise_sum_c2 = [];
ori_wise_med_c1 = []; ori_wise_sd_c1 = [];
ori_wise_med_c1_H = []; ori_wise_c1_H = [];
img_H = [];
fprintf('Starting HMAX\n')
count = 0;
for ver = 1:length(versions)
    verdir = [image_dir_base + versions(ver) + '/'];
    % loop across categories
    for cat = 1:5
        catdir = [verdir + cat_dict(5*(ver-1)+cat) + '/'];
        % get list of images in the folder of interest
        piclist = dir(catdir + '*.jpg');
        npics = length(piclist);
        
        % loop across pictures in folder, load and make grayscale
        for ipic = 1:npics
            % store version and category for each stimulus 
            stim_versions = [stim_versions; ver];
            stim_categories = [stim_categories; cat];
            
            % get the unique stimulus number that is used later
            stim_code = 1000 + cat_nums(5*(ver-1)+cat)*100+ipic;
            stim_codes = [stim_codes; stim_code];
            
            % load image
            pic_path = [catdir + piclist(ipic).name];
            pic_paths = [pic_paths; pic_path];
            picdat = load_pic(pic_path);
            
            % Extract layer responses
            [c2,c1] = extractC2forCell...
                (filters,filterSizes,c1Space,c1Scale,c1OL,patches,{picdat},nPatchSizes,patchSizes(1:3,:));
            
            % For each patch calculate responses
            patch_estim = get_patch_wise_estim(c1,c2);
            clear c2 
            % store patch-wise C1 and C2 values
            patch_wise_med_c1 = [patch_wise_med_c1; patch_estim.c1_median];
            patch_wise_med_c2 = [patch_wise_med_c2; patch_estim.c2_median];
            
            patch_wise_sd_c1 = [patch_wise_sd_c1; patch_estim.c1_sd];
            patch_wise_sd_c2 = [patch_wise_sd_c2; patch_estim.c2_sd];
            
            patch_wise_sum_c1 = [patch_wise_sum_c1; patch_estim.c1_sum];
            patch_wise_sum_c2 = [patch_wise_sum_c2; patch_estim.c2_sum];
            
            %c1_mean = [c1_mean; mean(patch_estim.c1_median)];
            %c2_mean = [c2_mean; mean(cell2mat(patch_estim.c2_median))];
            
            clear c1
            count = count + 1; 
            fprintf('\nFinished image %d. \n', count)
        end
    end
    
end

% save info
info = table(NaN(590,1), NaN(590,1), NaN(590,1), ... 
    'VariableNames', {'stimulus_code', 'version', 'category'});
info.stimulus_code = stim_codes;
info.version = stim_versions;
info.category = stim_categories; 
writetable(info, [save_dir_base + 'info.csv']) % as csv
save([save_dir_base + 'info.mat'], 'info')

% save all data
save([save_dir_base + 'HMAX_estims_allstims.mat'], 'info',...
    'patch_wise_med_c1', 'patch_wise_med_c2',...
    'patch_wise_sd_c1', 'patch_wise_sd_c2',...
    'patch_wise_sum_c1', 'patch_wise_sum_c2');



