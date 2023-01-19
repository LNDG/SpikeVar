clear
clc
% check if we are in the correct directory, change it if needed
cd('/Users/waschke/OneDrive/Old data and code/SV/OurCode/')
% add essential scripts
addpath('/Users/waschke/Documents/MATLAB/Code/Scripts/Essentials/')
addpath(genpath('/Users/waschke/Documents/Matlabtoolboxes/PLS/'))
% only needed for plotting
addpath(genpath('/Users/waschke/Documents/Matlabtoolboxes/RainCloudPlots/'))
addpath('/Users/waschke/Documents/Matlabtoolboxes/cbrewer/')
addpath('/Users/waschke/Documents/Matlabtoolboxes/Robust_Statistical_Toolbox/')

% VGG PCs
vggdir = '/Users/waschke/OneDrive/Old data and code/VGG16_output/';
load([vggdir 'VGG16_all_encoding_recog_layer1_sum_pc1.mat']);
load([vggdir 'VGG16_all_encoding_recog_layer2_sum_pc1.mat']);
load([vggdir 'VGG16_all_encoding_recog_layer3_sum_pc1.mat']);
load([vggdir 'VGG16_all_encoding_recog_layer4_sum_pc1.mat']);
load([vggdir 'VGG16_all_encoding_recog_layer5_sum_pc1.mat']);

load([vggdir 'VGG16_all_encoding_recog_layer1_wz_sd_pc1.mat']);
load([vggdir 'VGG16_all_encoding_recog_layer2_wz_sd_pc1.mat']);
load([vggdir 'VGG16_all_encoding_recog_layer3_wz_sd_pc1.mat']);
load([vggdir 'VGG16_all_encoding_recog_layer4_wz_sd_pc1.mat']);
load([vggdir 'VGG16_all_encoding_recog_layer5_wz_sd_pc1.mat']);

% also load HMAX components
hmax_dir = '/Users/waschke/OneDrive/Old data and code/HMAX output/';
load([hmax_dir 'HMAX_all_encoding_recog_C1_pc1.mat'])
load([hmax_dir 'HMAX_all_encoding_recog_C2_pc1.mat'])

% load spiking PE
load('SpikeVar_all_trials_incl_10_ms_PE.mat')
pe_table = perm_active_table;
% set of relevant IDs
cur_pat_ids = unique(pe_table.Participant);

%% loop across subjects
clear all_sub_neu_area
for sub_id = 1:length(cur_pat_ids)
    cid = cur_pat_ids(sub_id);
    % get relevant indices
    % set of neurons (all, hippocampus <3  or amygdala >2)
    % 1 = right Hippocampius, 2 = left Hippocampus
    % 3 = right Amygdala, 4 = left Amygdala
    num_tr_nums = [];
    for ises = 1:2
        % select neurons in amygdala (brainArea > 2)
        rel_ids = find(pe_table.Participant==cid &...
             pe_table.brainArea>2 &pe_table.Session==ises);

         % only continue if we have data
        if length(rel_ids)>=1
            % get number of trials (varies between participants and sessions)
            num_tr_nums(ises) = unique(pe_table.NumTrials(pe_table.Participant==cid &...
                 pe_table.brainArea>2 &pe_table.Session==ises));
            % save neuron info
            sub_neus{ises} = unique(pe_table.NeuronL(rel_ids));
            all_sub_neu_ids{sub_id, ises} = unique(pe_table.NeuronL(rel_ids));
            for ineu = 1:length(all_sub_neu_ids{sub_id, ises})
                all_sub_neu_area{sub_id, ises}(ineu) = unique(pe_table.brainArea(...
                    pe_table.NeuronL ==all_sub_neu_ids{sub_id, ises}(ineu)));
            end
            % extract PE for different motif lengths
            pe_dat_10_2{ises} = pe_table.psth10_permEN_2(rel_ids);
            pe_dat_10_3{ises} = pe_table.psth10_permEN_3(rel_ids);
            pe_dat_10_4{ises} = pe_table.psth10_permEN_4(rel_ids);
            % re-arrange data PE data
            re_mat_10_2{ises} = reshape(pe_dat_10_2{ises} , num_tr_nums(ises),...
                length(sub_neus{ises}));
            re_mat_10_3{ises} = reshape(pe_dat_10_3{ises} , num_tr_nums(ises),...
                length(sub_neus{ises}));
            re_mat_10_4{ises} = reshape(pe_dat_10_4{ises} , num_tr_nums(ises),...
                length(sub_neus{ises}));

            % kick out neurons that are very silent
            all_sub_neu_area{sub_id, ises}(:,(mean(re_mat_10_2{ises}))<.005) = [];
            % replicate neuron ids here so it can be used later
            all_sub_neu_area{sub_id, ises} = repmat(all_sub_neu_area{sub_id, ises},1,3);
            re_mat_10_2{ises}(:,(mean(re_mat_10_2{ises}))<.0001) = [];
            re_mat_10_3{ises}(:,(mean(re_mat_10_3{ises}))<.0001) = [];
            re_mat_10_4{ises}(:,(mean(re_mat_10_4{ises}))<.0001) = [];
            % save the final number of neurons per subject
            final_neu_num(sub_id,ises) = size(re_mat_10_2{ises},2);
            % also save average PE of remaining neurons
            avg_pe_mat_10_2{sub_id,ises} = mean(re_mat_10_2{ises});
            avg_pe_mat_10_3{sub_id,ises} = mean(re_mat_10_3{ises});
            avg_pe_mat_10_4{sub_id,ises} = mean(re_mat_10_4{ises});
            
            
            % another check if we still have neurons left
            if size(re_mat_10_2{ises},2)>1
                
                % get stimulus info for all layers
                % First, get PCs for the spatial sum of VGGG feature maps
                vgg1_inf = VGG1_sum_pc1_table.lay1_score(rel_ids,:);
                vgg2_inf = VGG2_sum_pc1_table.lay2_score(rel_ids,:);
                vgg3_inf = VGG3_sum_pc1_table.lay3_score(rel_ids,:);
                vgg4_inf = VGG4_sum_pc1_table.lay4_score(rel_ids,:);
                vgg5_inf = VGG5_sum_pc1_table.lay5_score(rel_ids,:);
                % next, the spatial SD of VGG feature maps
                vgg1_sd_inf = VGG1_sd_pc1_table.lay1_score(rel_ids,:);
                vgg2_sd_inf = VGG2_sd_pc1_table.lay2_score(rel_ids,:);
                vgg3_sd_inf = VGG3_sd_pc1_table.lay3_score(rel_ids,:);
                vgg4_sd_inf = VGG4_sd_pc1_table.lay4_score(rel_ids,:);
                vgg5_sd_inf = VGG5_sd_pc1_table.lay5_score(rel_ids,:);
                %  Next, for HMAX layers
                c1_inf = C1_pc1_table.C1_score(rel_ids);
                c2_inf = C2_pc1_table.C2_score(rel_ids);
                c1_sd_inf = C1_pc1_table.C1_sd_score(rel_ids);
                c2_sd_inf = C2_pc1_table.C2_sd_score(rel_ids);
                
                % re-arrange variables
                vgg1_inf_mat = reshape(vgg1_inf, num_tr_nums(ises),...
                    length(sub_neus{ises}), size(vgg1_inf,2));
                vgg2_inf_mat = reshape(vgg2_inf, num_tr_nums(ises),...
                    length(sub_neus{ises}), size(vgg2_inf,2));
                vgg3_inf_mat = reshape(vgg3_inf, num_tr_nums(ises),...
                    length(sub_neus{ises}), size(vgg3_inf,2));
                vgg4_inf_mat = reshape(vgg4_inf, num_tr_nums(ises),...
                    length(sub_neus{ises}), size(vgg4_inf,2));
                vgg5_inf_mat = reshape(vgg5_inf, num_tr_nums(ises),...
                    length(sub_neus{ises}), size(vgg5_inf,2));
                
                vgg1_sd_inf_mat = reshape(vgg1_sd_inf, num_tr_nums(ises),...
                    length(sub_neus{ises}), size(vgg1_sd_inf,2));
                vgg2_sd_inf_mat = reshape(vgg2_sd_inf, num_tr_nums(ises),...
                    length(sub_neus{ises}), size(vgg2_sd_inf,2));
                vgg3_sd_inf_mat = reshape(vgg3_sd_inf, num_tr_nums(ises),...
                    length(sub_neus{ises}), size(vgg3_sd_inf,2));
                vgg4_sd_inf_mat = reshape(vgg4_sd_inf, num_tr_nums(ises),...
                    length(sub_neus{ises}), size(vgg4_sd_inf,2));
                vgg5_sd_inf_mat = reshape(vgg5_sd_inf, num_tr_nums(ises),...
                    length(sub_neus{ises}), size(vgg5_sd_inf,2));
                
                c1_inf_mat = reshape(c1_inf, num_tr_nums(ises),...
                    length(sub_neus{ises}), size(c1_inf,2));
                c2_inf_mat = reshape(c2_inf, num_tr_nums(ises),...
                    length(sub_neus{ises}), size(c2_inf,2));
                c1_sd_inf_mat = reshape(c1_sd_inf, num_tr_nums(ises),...
                    length(sub_neus{ises}), size(c1_sd_inf,2));
                c2_sd_inf_mat = reshape(c2_sd_inf, num_tr_nums(ises),...
                    length(sub_neus{ises}), size(c2_sd_inf,2));
                
                % we'll only need one column of this
                vgg1_inf_rel{ises} = squeeze(vgg1_inf_mat(:,1,:));
                vgg2_inf_rel{ises} = squeeze(vgg2_inf_mat(:,1,:));
                vgg3_inf_rel{ises} = squeeze(vgg3_inf_mat(:,1,:));
                vgg4_inf_rel{ises} = squeeze(vgg4_inf_mat(:,1,:));
                vgg5_inf_rel{ises} = squeeze(vgg5_inf_mat(:,1,:));
                
                vgg1_sd_inf_rel{ises} = squeeze(vgg1_sd_inf_mat(:,1,:));
                vgg2_sd_inf_rel{ises} = squeeze(vgg2_sd_inf_mat(:,1,:));
                vgg3_sd_inf_rel{ises} = squeeze(vgg3_sd_inf_mat(:,1,:));
                vgg4_sd_inf_rel{ises} = squeeze(vgg4_sd_inf_mat(:,1,:));
                vgg5_sd_inf_rel{ises} = squeeze(vgg5_sd_inf_mat(:,1,:));
                
                c1_inf_rel{ises}   = squeeze(c1_inf_mat(:,1,:));
                c2_inf_rel{ises}   = squeeze(c2_inf_mat(:,1,:));
                c1_sd_inf_rel{ises}   = squeeze(c1_sd_inf_mat(:,1,:));
                c2_sd_inf_rel{ises}   = squeeze(c2_sd_inf_mat(:,1,:));
                

                % identify trials in which less than 1/3 of neurons are
                % spiking at all
                kick_ids{ises} = find(sum(sign(re_mat_10_2{ises}),2)<=...
                    (size(re_mat_10_2{ises},2)/3));
                
                % kick out these silent trials
                vgg1_inf_rel{ises}(kick_ids{ises},:) = [];
                vgg2_inf_rel{ises}(kick_ids{ises},:) = [];
                vgg3_inf_rel{ises}(kick_ids{ises},:) = [];
                vgg4_inf_rel{ises}(kick_ids{ises},:) = [];
                vgg5_inf_rel{ises}(kick_ids{ises},:) = [];
                
                vgg1_sd_inf_rel{ises}(kick_ids{ises},:) = [];
                vgg2_sd_inf_rel{ises}(kick_ids{ises},:) = [];
                vgg3_sd_inf_rel{ises}(kick_ids{ises},:) = [];
                vgg4_sd_inf_rel{ises}(kick_ids{ises},:) = [];
                vgg5_sd_inf_rel{ises}(kick_ids{ises},:) = [];
                
                c1_inf_rel{ises}(kick_ids{ises}) = [];
                c2_inf_rel{ises}(kick_ids{ises}) = [];
                c1_sd_inf_rel{ises}(kick_ids{ises}) = [];
                c2_sd_inf_rel{ises}(kick_ids{ises}) = [];
                
                % run PLS models
                % neural matrix consists of neuron and trial-wise PE for
                % motif lengths 2, 3, & 4
                datamat_lst = {[cell2mat(re_mat_10_2(ises)),...
                    cell2mat(re_mat_10_3(ises)),...
                    cell2mat(re_mat_10_4(ises))]};
                % set up some more options for PLS
                option = [];
                option.method = 3; % behavioural PLS
                option.num_perm = 1000; % number of permutations
                option.num_boot = 1000; % number of bootstraps
                option.cormode = 8; % rank correlation
                
                %  set up the stimulus info matrix for 
                % a model that includes all vgg and HMAX layers
                option.stacked_behavdata = [vgg1_inf_rel{ises}(:,1),...
                    vgg2_inf_rel{ises}(:,1),...
                    vgg3_inf_rel{ises}(:,1),...
                    vgg4_inf_rel{ises}(:,1)...
                    vgg5_inf_rel{ises}(:,1),...
                    vgg1_sd_inf_rel{ises}(:,1),...
                    vgg2_sd_inf_rel{ises}(:,1),...
                    vgg3_sd_inf_rel{ises}(:,1) vgg4_sd_inf_rel{ises}(:,1)...
                    vgg5_sd_inf_rel{ises}(:,1),...
                    c1_inf_rel{ises}, c2_inf_rel{ises}...
                    c1_sd_inf_rel{ises}, c2_sd_inf_rel{ises}]; %( 2-D numerical matrix )

                % kick out silent trials in the neural matrix
                datamat_lst{1}(kick_ids{ises},:) = [];
                % save neural data that has been used
                all_rel_mats{sub_id, ises} = datamat_lst{1};
                % save stimulus data that has been used
                all_rel_behav_dat{sub_id,ises} = option.stacked_behavdata;
                % number of trials
                num_trls_used = size(datamat_lst{1},1);
                % run model
                result_combined{ises} = pls_analysis(datamat_lst, num_trls_used, 1, option);
                
                % run one model using late layers
                option.stacked_behavdata = [vgg3_inf_rel{ises}(:,1),...
                    vgg4_inf_rel{ises}(:,1)...
                    vgg5_inf_rel{ises}(:,1),...
                    vgg3_sd_inf_rel{ises}(:,1) vgg4_sd_inf_rel{ises}(:,1)...
                    vgg5_sd_inf_rel{ises}(:,1),...
                    c2_inf_rel{ises}...
                    c2_sd_inf_rel{ises} ]; %( 2-D numerical matrix )
                result_late{ises} = pls_analysis(datamat_lst, num_trls_used, 1, option);
                
                % run one model using early layers
                option.stacked_behavdata = [vgg1_inf_rel{ises}(:,1),...
                    vgg2_inf_rel{ises}(:,1),...
                    vgg3_inf_rel{ises}(:,1),...
                    vgg1_sd_inf_rel{ises}(:,1),...
                    vgg2_sd_inf_rel{ises}(:,1),...
                    vgg3_sd_inf_rel{ises}(:,1),...
                    c1_inf_rel{ises}...
                    c1_sd_inf_rel{ises} ]; %( 2-D numerical matrix )
                result_early{ises} = pls_analysis(datamat_lst, num_trls_used, 1, option);
                

                % save subject-specific stimulus weights
                behav_weights.vggall{sub_id,ises} = result_combined{ises}.v;
                behav_weights.no12{sub_id,ises} = result_late{ises}.v;
                behav_weights.no45{sub_id,ises} = result_early{ises}.v;

                % save LVLV corrs
                lvlv_corrs.vggall{ises}(sub_id) = result_combined{ises}.lvlvcorr(1);
                lvlv_corrs.no12{ises}(sub_id) = result_late{ises}.lvlvcorr(1);
                lvlv_corrs.no45{ises}(sub_id) = result_early{ises}.lvlvcorr(1);

                
                % also get bounds of LVLV corrs  
                % use the actual lvlv corr distribution
                lvlv_corr_lims.vgall{ises}(sub_id,:) = [percentile(squeeze(...
                    result_combined{ises}.boot_result.lvlvcorrdistrib(1,:)),1),...
                    percentile(squeeze(...
                    result_combined{ises}.boot_result.lvlvcorrdistrib(1,:)),99)];
                lvlv_corr_lims.no12{ises}(sub_id,:) = [percentile(squeeze(...
                    result_late{ises}.boot_result.lvlvcorrdistrib(1,:)),1),...
                    percentile(squeeze(...
                    result_late{ises}.boot_result.lvlvcorrdistrib(1,:)),99)];
                lvlv_corr_lims.no45{ises}(sub_id,:) = [percentile(squeeze(...
                    result_early{ises}.boot_result.lvlvcorrdistrib(1,:)),1),...
                    percentile(squeeze(...
                    result_early{ises}.boot_result.lvlvcorrdistrib(1,:)),99)];

            end
        else
        end
    end
end
% save stuff
save('SpikeVar_early_vs_late_PE_PLS_amygdala_results.mat', 'behav_weights',...
    'lvlv_corrs','lvlv_corr_lims')
