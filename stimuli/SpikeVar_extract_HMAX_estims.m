clc
clear
% path to HMAX functions
addpath('/Users/kamp/PhD/spikevar/repo/SpikeVar/toolboxes/hmax-develop');
% paths to stimuli
% Use minimally renamed verion since it makes sorting things easier
image_dir_base = "/Users/kamp/PhD/spikevar/stimuli_rn";
save_dir_base = "/Users/kamp/PhD/spikevar/hmax_output";
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
all_HMAX_info = table(NaN(590,1),NaN(590,1),NaN(590,1), 'VariableNames', {'StimulusCode', 'C1mean', 'C2mean'});
stimcodes = [];
c1_vals = []; c2_vals = [];
patch_wise_c1 = []; patch_wise_c2 = [];
patch_wise_c1sd = []; patch_wise_c2sd = [];
ori_wise_med_c1 = []; ori_wise_sd_c1 = [];
ori_wise_med_c1_H = []; ori_wise_c1_H = [];
img_H = [];
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
            % load image, check if grayscale already
            temp_img = imread([catdir + piclist(ipic).name]);
            sz = size(temp_img);
            % either use gray scale image or transform RGB to grayscale
            % in either case, resize to 224*224
            % this is done to keep things comparable across images and
            % across methods since VGG16 only eats 224*224 images rn
            if size(sz)<3
                picdat{ipic} = double(imresize(temp_img,[224, 224]));
            else
                picdat{ipic} = double(imresize(rgb2gray(temp_img),[224,224]));
            end
            
            % For each patch calculate responses
            [c2,c1] = extractC2forCell...
                (filters,filterSizes,c1Space,c1Scale,c1OL,patches,   picdat(ipic),  nPatchSizes,patchSizes(1:3,:));
            c1_ori4 = cellfun(@(x) x(:,:,4), c1{1}, 'uni', 0); % take ori 4
            hmaxout.c1(:,:,ipic) = c1_ori4{8};
            hmaxout.c1median = cell2mat(cellfun(@(x) median(x(:)), c1_ori4, 'uni', 0)); % take median
            hmaxout.c2median  = cellfun(@median, c2, 'uni', 0); % take median
            
            hmaxout.c1SD = cell2mat(cellfun(@(x) std(x(:)), c1_ori4, 'uni', 0)); % take std
            hmaxout.c2SD  = cellfun(@std, c2, 'uni', 0); % take std
            
            hmaxout.picdat(:,:,ipic) = picdat{ipic};
            clear c2 
            % get the unique stimulus number that is used later
            stim_code = 1000+ cat_nums(5*(ver-1)+cat)*100+ipic;
            % append
            stimcodes = [stimcodes; stim_code];
            c1_vals = [c1_vals; mean(hmaxout.c1median)];
            c2_vals = [c2_vals; mean(cell2mat(hmaxout.c2median))];
            
            % also store patch-wise C1 and C2 values
            patch_wise_c1 = [patch_wise_c1; hmaxout.c1median];
            patch_wise_c2 = [patch_wise_c2; cell2mat(hmaxout.c2median)];
            
            % also store SD
            patch_wise_c1sd = [patch_wise_c1sd; hmaxout.c1SD];
            patch_wise_c2sd = [patch_wise_c2sd; cell2mat(hmaxout.c2SD)];
            
            % get entropy and orientation-wise estimates
            for k=1:length(c1Scale)-1
                for j=1:length(orientations)
                        median_c1_oris(k, j)=median(nonzeros(c1{1}{1, k}(:, :, j)));
                        sd_c1_oris(k, j)=std(nonzeros(c1{1}{1, k}(:, :, j)));
                        median_entropy_c1_oris(k, j) = median(nonzeros(entropyfilt((c1{1}{1, k}(:, :, j)))));
                        entropy_c1_oris(k,j) = entropy(c1{1}{1, k}(:,:,j));
                end
            end
            % entropy of image as well
            entropy_img = entropy(picdat{ipic});
            
            % append stuff
            ori_wise_med_c1(size(patch_wise_c1sd,1),:,:) = median_c1_oris; 
            ori_wise_sd_c1(size(patch_wise_c1sd,1),:,:) =  sd_c1_oris;
            ori_wise_med_c1_H(size(patch_wise_c1sd,1),:,:) =  median_entropy_c1_oris; 
            ori_wise_c1_H(size(patch_wise_c1sd,1),:,:) =  entropy_c1_oris;
            img_H(size(patch_wise_c1sd,1)) =  entropy_img;
            clear c1
        end
    end
    
end
% store all data in one table
all_HMAX_info.StimulusCode = stimcodes;
all_HMAX_info.C1mean = c1_vals;
all_HMAX_info.C2mean = c2_vals;

save([save_dir_base + 'HMAX_estims_allstims_April21_22.mat'], 'all_HMAX_info',...
    'patch_wise_c1', 'patch_wise_c2', 'patch_wise_c1sd', 'patch_wise_c2sd',...
    'ori_wise_med_c1', 'ori_wise_sd_c1', 'ori_wise_med_c1_H', 'ori_wise_c1_H',...
    'img_H');