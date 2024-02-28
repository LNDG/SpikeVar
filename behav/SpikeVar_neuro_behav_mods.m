clear
clc
% check if we are in the correct directory, change it if needed
cd('/Users/waschke/OneDrive/Old data and code/SV/OurCode/')
% add essential scripts
addpath('/Users/waschke/Documents/MATLAB/Code/Scripts/Essentials/')
addpath(genpath('/Users/waschke/Documents/Matlabtoolboxes/EntropyHub'))
addpath(genpath('/Users/waschke/Documents/Matlabtoolboxes/PLS/'))
addpath(genpath('/Users/waschke/Documents/Matlabtoolboxes/RainCloudPlots/'))
addpath('/Users/waschke/Documents/Matlabtoolboxes/cbrewer/')
addpath('/Users/waschke/Documents/Matlabtoolboxes/Robust_Statistical_Toolbox/')
% load spiking PE
load('SpikeVar_all_trials_incl_10_ms_PE.mat')
pe_table = perm_active_table;
cur_pat_ids = unique(pe_table.Participant);
% get task variant codes
for isub = 1:length(cur_pat_ids)
    for ises = 1:2
        if ~isempty(find(pe_table.Participant==...
                cur_pat_ids(isub)& pe_table.Session== ises))
            rel_it = find(pe_table.Participant==...
                cur_pat_ids(isub)& pe_table.Session== ises);
            task_vars(isub,ises) = unique(pe_table.TaskVariant(rel_it));
            image_dur(isub,ises) = unique(pe_table.learnImageDur(rel_it))/1000 ;
            avg_pe(isub,ises) = mean(pe_table.psth10_permEN_4(rel_it));
            spike_rate(isub,ises) =  mean(sum(pe_table.psth10(rel_it,1:image_dur(isub,ises)*100-1),2))/image_dur(isub,ises);
        else
            task_vars(isub,ises) = NaN;
            image_dur(isub,ises) = NaN;
            avg_pe(isub,ises) = NaN;
            spike_rate(isub,ises) = NaN;
        end
    end
end
% binarize task variant
task_vars(task_vars<3) = 0;
task_vars(task_vars ==3) = 1;
%% load PE results
load('SpikeVar_early_vs_late_PE_PLS_results')
all_stim_weights.PE = stim_weights;
all_lvlv_corrs.PE = lvlv_corrs;
all_lvlv_corr_lims.PE = lvlv_corr_lims;
%% load SPikenum results
load('SpikeVar_early_vs_late_Spikenum_PLS_results')
all_stim_weights.Spike = stim_weights;
all_lvlv_corrs.Spike = lvlv_corrs;
all_lvlv_corr_lims.Spike = lvlv_corr_lims;
%% load Amygdala results
load('SpikeVar_early_vs_late_PE_PLS_amygdala_results')
all_stim_weights.Amyg_PE = stim_weights;
all_lvlv_corrs.Amyg_PE = lvlv_corrs;
all_lvlv_corr_lims.Amyg_PE = lvlv_corr_lims;
%% load peformance info
load('SpikeVar_all_behav_info.mat')
%% Prepare linear models
% rank - regression due to small sample 
rel_inds = find(all_lvlv_corrs.PE.vggall{1}~=0);
s1_rank_tab = table(tiedrank(all_lvlv_corrs.PE.no45{1}(rel_inds))',...
    tiedrank(all_lvlv_corrs.PE.no12{1}(rel_inds))',...
    tiedrank(all_lvlv_corrs.PE.vggall{1}(rel_inds))',...
    tiedrank(all_lvlv_corrs.Spike.no12{1}(rel_inds))',...
    tiedrank(all_lvlv_corrs.Amyg_PE.no12{1}(rel_inds))',...
    tiedrank(all_lvlv_corrs.Amyg_PE.no45{1}(rel_inds))',...
    tiedrank(avg_behav_dat{1}.rlogP(rel_inds)),...
    tiedrank(avg_behav_dat{1}.rcwP(rel_inds)),...
    tiedrank(avg_behav_dat{1}.dprim(rel_inds)),...
    tiedrank(avg_behav_dat{1}.rConf(rel_inds)),...
    tiedrank(avg_behav_dat{1}.crit(rel_inds)),...
    tiedrank(all_neu_nums(rel_inds,1)),...
    tiedrank(all_trl_nums(rel_inds,1)),...
    (avg_behav_dat{1}.Age(rel_inds)),...
    categorical(task_vars(rel_inds,1)),...
    categorical(image_dur(rel_inds,1)),...
    tiedrank(avg_pe(rel_inds,1)),...
    tiedrank(spike_rate(rel_inds,1)),...
    'VariableNames', {'PE_n45','PE_n12','PE_A','Spike_n12','Amyg_n12',...
    'Amyg_n45','rlogP', 'rcwP', 'dprim', 'rconf',...
    'crit','numNeus', 'numTrls','Age', 'TaskV', 'Dur', 'Avg_PE', 'Spike_rate'});
%% extra table for amaygdala controls
% We only need this since the sample 
rel_inds = find(all_lvlv_corrs.PE.vggall{1}~=0 & all_lvlv_corrs.Amyg_PE.no12{1}~=0);
s1_rank_amyg_tab = table(tiedrank(all_lvlv_corrs.PE.no45{1}(rel_inds))',...
    tiedrank(all_lvlv_corrs.PE.no12{1}(rel_inds))',...
    tiedrank(all_lvlv_corrs.PE.vggall{1}(rel_inds))',...
    tiedrank(all_lvlv_corrs.Spike.no12{1}(rel_inds))',...
    tiedrank(all_lvlv_corrs.Amyg_PE.no12{1}(rel_inds))',...
    tiedrank(all_lvlv_corrs.Amyg_PE.no45{1}(rel_inds))',...
    tiedrank(avg_behav_dat{1}.rlogP(rel_inds)),...
    tiedrank(avg_behav_dat{1}.rcwP(rel_inds)),...
    tiedrank(avg_behav_dat{1}.dprim(rel_inds)),...
    tiedrank(avg_behav_dat{1}.rConf(rel_inds)),...
    tiedrank(avg_behav_dat{1}.crit(rel_inds)),...
    tiedrank(all_neu_nums(rel_inds,1)),...
    tiedrank(all_trl_nums(rel_inds,1)),...
    (avg_behav_dat{1}.Age(rel_inds)),...
    categorical(task_vars(rel_inds,1)),...
    categorical(image_dur(rel_inds,1)),...
    tiedrank(avg_pe(rel_inds,1)),...
    tiedrank(spike_rate(rel_inds,1)),...
    'VariableNames', {'PE_n45','PE_n12','PE_A','Spike_n12','Amyg_n12',...
    'Amyg_n45','rlogP', 'rcwP', 'dprim', 'rconf',...
    'crit','numNeus', 'numTrls','Age', 'TaskV', 'Dur', 'Avg_PE', 'Spike_rate'});

%% PCA on performance metrics
behav_metrics = table2array(s1_rank_tab(:,7:10));
[behav_coeff,behav_score,...
    behav_latent,~,behav_exp_var,~] =...
    pca(zscore(behav_metrics), 'Centered', true,'VariableWeights', 'Variance');
loadings = (corr(behav_score,...
    behav_metrics(:,1:size(behav_score,2))))';
exp_var_1st_pc = behav_exp_var(1);
exp_normvar_1st_pc = behav_exp_var(1)/...
    (100/size(behav_exp_var,1));

spss_scores = [-.49096 -1.03220 -.95593 -1.34061 -.99076 1.58337 .08288,...
    -.88392 -1.54987 .25994 -.72726 .07078 -1.24414 .23109 .13330 .96957,...
    1.60498 -1.65079 -.11430 .85260 .91304 .15395 -1.07072 -.75862 1.30807,...
    .69051 .60916 -.16001 -.79968 1.65848 .52723 .34587 1.77494];
corr(behav_score(:,1),spss_scores');

s1_rank_tab.behav_pc = behav_score(:,1);
s1_rank_tab.rank_behav_pc = tiedrank(behav_score(:,1));

rel_inds_hipp = find(all_lvlv_corrs.PE.no12{1}~=0);
rel_inds_amyg = find(all_lvlv_corrs.Amyg_PE.no12{1}(rel_inds_hipp)~=0);
s1_rank_amyg_tab.behav_pc = behav_score(rel_inds_amyg,1);
s1_rank_amyg_tab.rank_behav_pc = tiedrank(behav_score(rel_inds_amyg,1));
%% Linear models: each layer against performance PC
for imod = 1:3
    PE_perf_model_ranks{imod} = fitlm(s1_rank_tab, ['rank_behav_pc ~  ' ...
        s1_rank_tab.Properties.VariableNames{imod}]);
    
    PE_mod_eval.estims.behav_pc(imod) = PE_perf_model_ranks{imod}.Coefficients.Estimate(2);
    
    PE_mod_eval.tval.behav_pc(imod) = PE_perf_model_ranks{imod}.Coefficients.tStat(2);
    
    PE_mod_eval.r2.behav_pc(imod) = PE_perf_model_ranks{imod}.Rsquared.Adjusted;
    
    PE_mod_eval.p.behav_pc(imod) = PE_perf_model_ranks{imod}.Coefficients.pValue(2);
    
    temp = coefCI(PE_perf_model_ranks{imod}); PE_mod_eval.CI.behav_pc(imod,:) = temp(2,:);
end
%% plot
cols.pe = [52, 160, 164;...
    26, 117, 159;...
    24, 78, 119]/255;
cols.spikenum = [220, 47, 2;...
    244, 140, 6;...
    255, 186, 8]/255;

figure
xvals= [1:3];

for b = 1:size(PE_mod_eval.estims.behav_pc,2)
    b_part_pe(b) = bar(xvals(b),PE_mod_eval.estims.behav_pc(b));
    hold on
    plot(repmat(xvals(b),2,1), [PE_mod_eval.CI.behav_pc(b,1),...
        PE_mod_eval.CI.behav_pc(b,2)]', ...
        'Color',[.7 .7 .7] );
    b_part_pe(b).FaceColor = cols.pe(b,:);
    b_part_pe(b).BarWidth = .5;
    b_part_pe(b).FaceAlpha = 1;
    
    
end

ylabel('partial eta')
xticks([1:3])
xticklabels({'Early layers',  'Late layers', 'All layers'})

labs = {'Early layers',  'Late layers', 'All layers'};
% scatters
figure
for imod = 1:3
    subplot(3,1,imod)
    scatter(table2array(s1_rank_tab(:,imod)), s1_rank_tab.rank_behav_pc,140,...
        'filled', 'MarkerFaceAlpha', .8, 'MarkerFaceColor', cols.pe(imod,:))
    lsline
    xticks([])
    xlabel('PE modulation [ranks]')
    title(labs{imod})
    ylabel('Performance [ranks]')
    xxs = xlim(); yys = ylim();
    text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2),['partial eta = ' num2str(...
        round(PE_mod_eval.estims.behav_pc(imod)*1000)/1000)...
        ' p = ' num2str(round(PE_mod_eval.p.behav_pc(imod)*1000)/1000)])
end
%% Contrast layer-specific explanatory power
% model performance as a function of all layers first
PE_perf_model_early_late_ranks = fitlm(s1_rank_tab, ['rank_behav_pc ~  ' ...
    'PE_n12 + PE_n45 ']);
PE_perf_model_all_late_ranks = fitlm(s1_rank_tab, ['rank_behav_pc ~  ' ...
    'PE_n12 + PE_A ']);
% --> Late layers survive as only predictor in both cases

% get partial etas
PE_perf_cont_early_mod_ = fitlm(s1_rank_tab, ['rank_behav_pc ~  ' ...
    'PE_n45 ']);
PE_perf_cont_all_mod_ = fitlm(s1_rank_tab, ['rank_behav_pc ~  ' ...
    'PE_A ']);
PE_pe_cont_early_mod_ = fitlm(s1_rank_tab, ['PE_n12 ~  ' ...
    'PE_n45 ']);
PE_pe_cont_all_mod_ = fitlm(s1_rank_tab, ['PE_n12 ~  ' ...
    'PE_A ']);
PE_compare.eta.early_vs_late = (corr(tiedrank(PE_perf_cont_early_mod_.Residuals.Raw),...
    tiedrank(PE_pe_cont_early_mod_.Residuals.Raw)));
PE_compare.eta.all_vs_late = (corr(tiedrank(PE_perf_cont_all_mod_.Residuals.Raw),...
    tiedrank(PE_pe_cont_all_mod_.Residuals.Raw)));
% get correct p-value
% first transform r to t
% we need the "real" n
n_actual = size(s1_rank_tab,1)-...
    (size(PE_perf_model_early_late_ranks.Coefficients,1)-2);
PE_compare.t.early_late = r_to_t(PE_compare.eta.early_vs_late, n_actual);
PE_compare.p.early_late = t_to_p(PE_compare.t.early_late, n_actual);

PE_compare.t.all_late = r_to_t(PE_compare.eta.all_vs_late, n_actual);
PE_compare.p.all_late = t_to_p(PE_compare.t.all_late, n_actual);

% Boostrap CI
mycorr = @(x1,x2) corr(x1,x2,'Type','Pearson');
nI = 10000;
[PE_compare.etaCI.early_vs_late , ~] = bootci(nI,{mycorr,tiedrank(PE_perf_cont_early_mod_.Residuals.Raw),...
    tiedrank(PE_pe_cont_early_mod_.Residuals.Raw)});
[PE_compare.etaCI.all_vs_late , ~] = bootci(nI,{mycorr,tiedrank(PE_perf_cont_all_mod_.Residuals.Raw),...
    tiedrank(PE_pe_cont_all_mod_.Residuals.Raw)});
%% control late layer PE with late layer Spike Number
PE_perf_model_spike_late_ranks = fitlm(s1_rank_tab, ['rank_behav_pc ~  ' ...
    'PE_n12 + Spike_n12 ']);
% --> late layer PE survives

% Get partial etas
PE_perf_cont_spiken_mod = fitlm(s1_rank_tab, ['rank_behav_pc ~  ' ...
    'Spike_n12 ']);
PE_pe_cont_spiken_mod = fitlm(s1_rank_tab, ['PE_n12 ~  ' ...
    'Spike_n12 ']);

PE_compare.eta.spikenum_vs_PElate = (corr(tiedrank(PE_perf_cont_spiken_mod.Residuals.Raw),...
    tiedrank(PE_pe_cont_spiken_mod.Residuals.Raw)));
% get correct p-value
% first transform r to t
% we need the "real" n
n_actual = size(s1_rank_tab,1)-...
    (size(PE_perf_model_spike_late_ranks.Coefficients,1)-2);
PE_compare.t.spikenum_vs_PElate = r_to_t(PE_compare.eta.spikenum_vs_PElate, n_actual);
PE_compare.p.spikenum_vs_PElate = t_to_p(PE_compare.t.spikenum_vs_PElate, n_actual);

% Boostrap CI
mycorr = @(x1,x2) corr(x1,x2,'Type','Pearson');
nI = 10000;
[PE_compare.etaCI.spikenum_vs_PElate , ~] = bootci(nI,{mycorr,tiedrank(PE_perf_cont_spiken_mod.Residuals.Raw),...
    tiedrank(PE_pe_cont_spiken_mod.Residuals.Raw)});

%% control late layer PE for amygdala PE
% find participants without amygdala neurons

PE_perf_model_amyg_late_ranks = fitlm(s1_rank_amyg_tab, ['rank_behav_pc ~  ' ...
    'PE_n12 + Amyg_n12 + Amyg_n45']);
% --> late layer PE survives

% Get partial etas
PE_perf_cont_amyg_mod = fitlm(s1_rank_amyg_tab, ['rank_behav_pc ~  ' ...
    'Amyg_n12 + Amyg_n45']);
PE_pe_cont_amyg_mod = fitlm(s1_rank_amyg_tab, ['PE_n12 ~  ' ...
    'Amyg_n12 + Amyg_n45']);

PE_compare.eta.AmygLate_vs_PElate = (corr(tiedrank(PE_perf_cont_amyg_mod.Residuals.Raw),...
    tiedrank(PE_pe_cont_amyg_mod.Residuals.Raw)));
% get correct p-value
% first transform r to t
% we need the "real" n
n_actual = size(s1_rank_amyg_tab,1)-...
    (size(PE_perf_model_amyg_late_ranks.Coefficients,1)-2);
PE_compare.t.AmygLate_vs_PElate = r_to_t(PE_compare.eta.AmygLate_vs_PElate, n_actual);
PE_compare.p.AmygLate_vs_PElate = t_to_p(PE_compare.t.AmygLate_vs_PElate, n_actual);

% Boostrap CI
mycorr = @(x1,x2) corr(x1,x2,'Type','Pearson');
nI = 10000;
[PE_compare.etaCI.AmygLate_vs_PElate , ~] = bootci(nI,{mycorr,tiedrank(PE_perf_cont_amyg_mod.Residuals.Raw),...
    tiedrank(PE_pe_cont_amyg_mod.Residuals.Raw)});

%% control late PE for additional potential confounders
PE_perf_model_late_controls_ranks = fitlm(s1_rank_tab, ['rank_behav_pc ~  ' ...
    'PE_n12 +  numNeus + numTrls+Age+TaskV+Dur']);
% --> PE survives

% Get partial etas
PE_perf_cont_control_mod = fitlm(s1_rank_tab, ['rank_behav_pc ~  ' ...
    'numNeus + numTrls+Age+TaskV+Dur']);
PE_pe_cont_control_mod_ = fitlm(s1_rank_tab, ['PE_n12 ~  ' ...
    'numNeus + numTrls+Age+TaskV+Dur ']);

PE_compare.eta.PElate_controls = (corr(tiedrank(PE_perf_cont_control_mod.Residuals.Raw),...
    tiedrank(PE_pe_cont_control_mod_.Residuals.Raw)));
% get correct p-value
% first transform r to t
% we need the "real" n
n_actual = size(s1_rank_tab,1)-...
    (size(PE_perf_cont_control_mod.Coefficients,1)-2);
PE_compare.t.PElate_controls = r_to_t(PE_compare.eta.PElate_controls, n_actual);
PE_compare.p.PElate_controls = t_to_p(PE_compare.t.PElate_controls, n_actual);

% Boostrap CI
mycorr = @(x1,x2) corr(x1,x2,'Type','Pearson');
nI = 10000;
[PE_compare.etaCI.PElate_controls , ~] = bootci(nI,{mycorr,tiedrank(PE_perf_cont_control_mod.Residuals.Raw),...
    tiedrank(PE_pe_cont_control_mod_.Residuals.Raw)});
%% Plots
cols.pe = [52, 160, 164;...
    26, 117, 159;...
    24, 78, 119]/255;
figure
subplot(2,2,1)
hold on
scatter(table2array(s1_rank_tab(:,2)), s1_rank_tab.rank_behav_pc,140,...
    'filled', 'MarkerFaceAlpha', .8, 'MarkerFaceColor', cols.pe(2,:))
lsline
xticks([])
xlabel('PE modulation [ranks]')
title('Zero order')
ylabel('Performance [ranks]')
xxs = xlim(); yys = ylim();
text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2),['partial eta = ' num2str(...
    round(PE_mod_eval.estims.behav_pc(2)*1000)/1000)...
    ' p = ' num2str(round(PE_mod_eval.p.behav_pc(2)*1000)/1000)])

subplot(2,2,2)
scatter(tiedrank(PE_perf_cont_early_mod_.Residuals.Raw),...
    tiedrank(PE_pe_cont_early_mod_.Residuals.Raw),140,...
    'filled', 'MarkerFaceAlpha', 1, 'MarkerFaceColor', cols.pe(2,:))
lsline
xticks([])
xlabel('PE modulation [ranks]')
title('Controlled for early layers')
ylabel('Performance [ranks]')
xxs = xlim(); yys = ylim();
text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2),['partial eta = ' num2str(...
    round(PE_compare.eta.early_vs_late*1000)/1000)...
    ' p = ' num2str(round(PE_compare.p.early_late*1000)/1000)])

subplot(2,2,3)
scatter(tiedrank(PE_perf_cont_spiken_mod.Residuals.Raw),...
    tiedrank(PE_pe_cont_spiken_mod.Residuals.Raw),140,...
    'filled', 'MarkerFaceAlpha', 1, 'MarkerFaceColor', cols.pe(2,:))
lsline
xticks([])
xlabel('PE modulation [ranks]')
title('Controlled for SpikeNum')
ylabel('Performance [ranks]')
xxs = xlim(); yys = ylim();
text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2),['partial eta = ' num2str(...
    round(PE_compare.eta.spikenum_vs_PElate*1000)/1000)...
    ' p = ' num2str(round(PE_compare.p.spikenum_vs_PElate*10000)/10000)])
subplot(2,2,4)
scatter(tiedrank(PE_perf_cont_control_mod.Residuals.Raw),...
    tiedrank(PE_pe_cont_control_mod_.Residuals.Raw),140,...
    'filled', 'MarkerFaceAlpha', 1, 'MarkerFaceColor', cols.pe(2,:))
lsline
xticks([])
xlabel('PE modulation [ranks]')
title('Controlled for potential confounders')
ylabel('Performance [ranks]')
xxs = xlim(); yys = ylim();
text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2),['partial eta = ' num2str(...
    round(PE_compare.eta.PElate_controls*1000)/1000)...
    ' p = ' num2str(round(PE_compare.p.PElate_controls*10000)/10000)])

%long horizontal line of bar graphs

pl_etas = [PE_mod_eval.estims.behav_pc(1),PE_mod_eval.estims.behav_pc(2),...
    PE_mod_eval.estims.behav_pc(3),...
    PE_compare.eta.early_vs_late, PE_compare.eta.all_vs_late,...
    PE_compare.eta.spikenum_vs_PElate,...
    PE_compare.eta.PElate_controls];
pl_CIs = [PE_mod_eval.CI.behav_pc(1,:)',PE_mod_eval.CI.behav_pc(2,:)'...
    PE_mod_eval.CI.behav_pc(3,:)',...
    PE_compare.etaCI.early_vs_late, PE_compare.etaCI.all_vs_late,...
    PE_compare.etaCI.spikenum_vs_PElate,...
    PE_compare.etaCI.PElate_controls];


xvals= [.5, 1, 1.5,...
    2.5,3,...
    4,...
    5];
figure
for b = 1:length(pl_etas)
    b_part_pe(b) = bar(xvals(b),pl_etas(b));
    hold on
    plot(repmat(xvals(b),2,1), [pl_CIs(1,b),...
        pl_CIs(2,b)]', ...
        'Color',[.8 .8 .8] );
    b_part_pe(b).FaceColor = cols.pe(2,:);
    b_part_pe(b).BarWidth = .5;
    if b >1
        b_part_pe(b).FaceAlpha = 1/2;
    end
end
xticks(xvals)
xticklabels({'Early', 'Late', 'All',...
    'Cont early', 'Cont all',...
    'Cont SpikeNum',...
    'Cont Confounds'})
ylabel('Partial eta')

%% zero order scatters for amygala
zero_order_amyg_mod = fitlm(s1_rank_amyg_tab, ['rank_behav_pc ~  ' ...
    'Amyg_n45 ']);
figure
hold on
scatter(s1_rank_amyg_tab.Amyg_n12, s1_rank_amyg_tab.rank_behav_pc,140,...
    'filled', 'MarkerFaceAlpha', .8, 'MarkerFaceColor', cols.pe(2,:))
lsline
xticks([])
xlabel('Amygdala PE modulation [ranks]')
title('Zero order')
ylabel('Performance [ranks]')
xxs = xlim(); yys = ylim();
text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2),['partial eta = ' num2str(...
    round(zero_order_amyg_mod.Coefficients.Estimate(2)*1000)/1000)...
    ' p = ' num2str(round(zero_order_amyg_mod.Coefficients.pValue(2)*1000)/1000)])
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ANYTHING FROM HERE ON ESSENTIALLY IS SUPPLEMENTAL MATERIAL %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% linear models
for imod = 1:3
    PE_dprime_model_ranks{imod} = fitlm(s1_rank_tab, ['dprim ~  ' ...
        s1_rank_tab.Properties.VariableNames{imod}]);
    PE_logp_model_ranks{imod} = fitlm(s1_rank_tab, ['rlogP ~  ' ...
        s1_rank_tab.Properties.VariableNames{imod}]);
    PE_rcwp_model_ranks{imod} = fitlm(s1_rank_tab, ['rcwP ~  ' ...
        s1_rank_tab.Properties.VariableNames{imod}]);
    PE_rconf_model_ranks{imod} = fitlm(s1_rank_tab, ['rconf ~  ' ...
        s1_rank_tab.Properties.VariableNames{imod}]);
    PE_crit_model_ranks{imod} = fitlm(s1_rank_tab, ['crit ~  ' ...
        s1_rank_tab.Properties.VariableNames{imod}]);
    
    PE_mod_eval.estims.dprime(imod) = PE_dprime_model_ranks{imod}.Coefficients.Estimate(2);
    PE_mod_eval.estims.rlogp(imod) = PE_logp_model_ranks{imod}.Coefficients.Estimate(2);
    PE_mod_eval.estims.rcwp(imod) = PE_rcwp_model_ranks{imod}.Coefficients.Estimate(2);
    PE_mod_eval.estims.rconf(imod) = PE_rconf_model_ranks{imod}.Coefficients.Estimate(2);
    PE_mod_eval.estims.crit(imod) = PE_crit_model_ranks{imod}.Coefficients.Estimate(2);
    
    PE_mod_eval.tval.dprime(imod) = PE_dprime_model_ranks{imod}.Coefficients.tStat(2);
    PE_mod_eval.tval.rlogp(imod) = PE_logp_model_ranks{imod}.Coefficients.tStat(2);
    PE_mod_eval.tval.rcwp(imod) = PE_rcwp_model_ranks{imod}.Coefficients.tStat(2);
    PE_mod_eval.tval.rconf(imod) = PE_rconf_model_ranks{imod}.Coefficients.tStat(2);
    PE_mod_eval.tval.crit(imod) = PE_crit_model_ranks{imod}.Coefficients.tStat(2);
    
    PE_mod_eval.r2.dprime(imod) = PE_dprime_model_ranks{imod}.Rsquared.Adjusted;
    PE_mod_eval.r2.rlogp(imod) = PE_logp_model_ranks{imod}.Rsquared.Adjusted;
    PE_mod_eval.r2.rcwp(imod) = PE_rcwp_model_ranks{imod}.Rsquared.Adjusted;
    PE_mod_eval.r2.rconf(imod) = PE_rconf_model_ranks{imod}.Rsquared.Adjusted;
    PE_mod_eval.r2.crit(imod) = PE_crit_model_ranks{imod}.Rsquared.Adjusted;
    
    PE_mod_eval.p.dprime(imod) = PE_dprime_model_ranks{imod}.Coefficients.pValue(2);
    PE_mod_eval.p.rlogp(imod) = PE_logp_model_ranks{imod}.Coefficients.pValue(2);
    PE_mod_eval.p.rcwp(imod) = PE_rcwp_model_ranks{imod}.Coefficients.pValue(2);
    PE_mod_eval.p.rconf(imod) = PE_rconf_model_ranks{imod}.Coefficients.pValue(2);
    PE_mod_eval.p.crit(imod) = PE_crit_model_ranks{imod}.Coefficients.pValue(2);
    
    temp = coefCI(PE_dprime_model_ranks{imod}); PE_mod_eval.CI.dprime(imod,:) = temp(2,:);
    temp = coefCI(PE_logp_model_ranks{imod}); PE_mod_eval.CI.rlogp(imod,:) = temp(2,:);
    temp = coefCI(PE_rcwp_model_ranks{imod});PE_mod_eval.CI.rcwp(imod,:) = temp(2,:);
    temp = coefCI(PE_rconf_model_ranks{imod});PE_mod_eval.CI.rconf(imod,:) = temp(2,:);
    temp = coefCI(PE_crit_model_ranks{imod});PE_mod_eval.CI.crit(imod,:) = temp(2,:);
end

%% compare early and late layer effects
early_late_dprime_mod = fitlm(s1_rank_tab, ['dprim ~  ' ...
    s1_rank_tab.Properties.VariableNames{1} '+' s1_rank_tab.Properties.VariableNames{2}]);
early_late_rlogp_mod = fitlm(s1_rank_tab, ['rlogP ~  ' ...
    s1_rank_tab.Properties.VariableNames{1} '+' s1_rank_tab.Properties.VariableNames{2}]);
early_late_rcwp_mod = fitlm(s1_rank_tab, ['rcwP ~  ' ...
    s1_rank_tab.Properties.VariableNames{1} '+' s1_rank_tab.Properties.VariableNames{2}]);
early_late_conf_mod = fitlm(s1_rank_tab, ['rconf ~  ' ...
    s1_rank_tab.Properties.VariableNames{1} '+' s1_rank_tab.Properties.VariableNames{2}]);
early_late_crit_mod = fitlm(s1_rank_tab, ['crit ~  ' ...
    s1_rank_tab.Properties.VariableNames{1} '+' s1_rank_tab.Properties.VariableNames{2}]);
%% Plot partial etas including CI
estim_fbar = [PE_mod_eval.estims.rlogp;...
    PE_mod_eval.estims.dprime;...
    PE_mod_eval.estims.rconf;...
    PE_mod_eval.estims.rcwp;...
    PE_mod_eval.estims.crit];

estim_ci_fbar.low = [PE_mod_eval.CI.rlogp(:,1)';...
    PE_mod_eval.CI.dprime(:,1)';...
    PE_mod_eval.CI.rconf(:,1)';...
    PE_mod_eval.CI.rcwp(:,1)';...
    PE_mod_eval.CI.crit(:,1)'];

estim_ci_fbar.hi = [PE_mod_eval.CI.rlogp(:,2)';...
    PE_mod_eval.CI.dprime(:,2)';...
    PE_mod_eval.CI.rconf(:,2)';...
    PE_mod_eval.CI.rcwp(:,2)';...
    PE_mod_eval.CI.crit(:,2)'];
figure
b_r2 = bar(estim_fbar);

ylabel('partial eta')
xticklabels({'Accuracy', 'dprime', 'Confidence', 'Conf-wghtd-Acc',...
    'Criterion'})

b_r2(1).FaceColor = [100, 58, 113]/255;
b_r2(2).FaceColor = [139, 95, 191]/255;
b_r2(3).FaceColor = [209, 131, 201]/255;
title('Zero order model results')

% add error bars
hold on
ngroups = size(estim_fbar, 1);
nbars = size(estim_fbar, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    plot(repmat(x,2,1), [estim_ci_fbar.low(:,i), estim_ci_fbar.hi(:,i)]', ...
        'Color',[.7 .7 .7] );
end
legend('Early layers', 'Late layers', 'All layers')
%% check for seizure side vs recording side effects
% NOsessions.diagnosisCode (0 to 8):
% a code indicating the location of epileptic focal point:
% 0) not localized,
%1) Right Mesial Temporal,
%2) Left Mesial Temporal,
%3) Right Neocortical Temporal,
%4) Left Neocortical Temporal,
%5) Right Lateral Frontal,
%6) Left Lateral Frontal,
%7) Bilateral Independent Temporal,
%8) Bilateral Independent Frontal,
%9) Right Other,
%10) Left Other.
% brainAreaOfCell:  
%1=Right Hippocampus, 
%2=Left Hippocampus, 
%3=Right Amygdala, 
%4=Left Amygdala

for isub = 1:length(cur_pat_ids)
    for ises = 1:2
        if ~isempty(find(pe_table.Participant==...
                cur_pat_ids(isub)& pe_table.Session== ises))
            rel_it = find(pe_table.Participant==...
                cur_pat_ids(isub)& pe_table.Session== ises);
            SOZ(isub,ises) = unique(pe_table.DiagnosisCode(rel_it));
            if ismember(SOZ(isub,ises), [1,3,7])
                SOZ_side(isub,ises) = 0;
            elseif ismember(SOZ(isub,ises), [2,4,7])
                SOZ_side(isub,ises) = 1;
            else 
                SOZ_side(isub,ises) = NaN;
            end
            % recording side
            temp_loc = unique(pe_table.brainArea(rel_it));
            hipp_loc = temp_loc(ismember(temp_loc, [1,2]));
            if length(hipp_loc) ==1
                % only one side: 0 = right, 1 = left
                rec_loc(isub,ises) = hipp_loc-1;
            else
                rec_loc(isub,ises) = 2;
            end
           
        else
            SOZ(isub,ises) = NaN;
            rec_loc(isub,ises) = NaN;
            SOZ_side(isub,ises) = NaN;
        end
        % check matches
        if SOZ_side(isub,ises) == 0 && ismember(rec_loc(isub,ises), [0,2])
            side_same(isub,ises)  = 1;
        elseif SOZ_side(isub,ises) == 1 && ismember(rec_loc(isub,ises), [1,2])
            side_same(isub,ises)  = 1;
        else 
            side_same(isub,ises)  = 0;
        end
    end
end
%%
s1_rank_tab.SOZ_inside = side_same(rel_inds,1);

PE_SOZ_con_mod = fitlm(s1_rank_tab, ['rank_behav_pc ~  ' ...
    'SOZ_inside']);
