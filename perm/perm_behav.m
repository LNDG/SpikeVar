clear 
clc
base_dir = "XXXX"; % enter the directory of the repository here
cd(base_dir)
data_dir = fullfile(base_dir , "data/");
neuro_dir = fullfile(base_dir, "output", "neuro/");
behav_dir = fullfile(base_dir, "output", "behav/");
perm_dir = fullfile(base_dir, "output", "perm/");
addpath(genpath(base_dir));

compute_estimates = 0; 
compute_permstats = 0;

%% load spike data and pls results
% spike data
load([neuro_dir + 'SpikeVar_spike_data_PE_or.mat'], 'spike_data_table')
spike_summary = c_get_spike_summay(spike_data_table, neuro_dir);

% pls results
load([perm_dir + 'pls_permutation_results_or.mat']);

% neuron and trial nums
neuron_nums = arrayfun(@(x) x.neuron_nums, info.hippocampus);
trial_nums = arrayfun(@(x) x.trial_nums, info.hippocampus);

% relevant indices
rel_inds.hpc = find(pls_results.pe.hippocampus.all_layers(:,1)~=0);
rel_inds.amy = find(pls_results.pe.hippocampus.all_layers(:,1)~=0 & ... 
    pls_results.pe.amygdala.all_layers(:,1)~=0);

%% load behav data & rank
load([data_dir + 'SpikeVar_learn_recog_info.mat'], 'fulltable_learn', 'fulltable_recog');
pat_ids = unique(spike_data_table.Participant);
behav_dat = a_extract_behav_info(fulltable_learn, fulltable_recog, pat_ids, behav_dir);

% hippocampus
% convert behav to rank and compute pc 
rank_table.hpc = table( ... 
    tiedrank(behav_dat.session_1.rlogP(rel_inds.hpc)),...
    tiedrank(behav_dat.session_1.rcwP(rel_inds.hpc)),...
    tiedrank(behav_dat.session_1.dprim(rel_inds.hpc)),...
    tiedrank(behav_dat.session_1.rConf(rel_inds.hpc)),...
    tiedrank(neuron_nums(rel_inds.hpc,1)),...
    tiedrank(trial_nums(rel_inds.hpc,1)),...
    (behav_dat.session_1.Age(rel_inds.hpc)),...
    categorical(spike_summary.task_vars(rel_inds.hpc,1)),...
    categorical(spike_summary.image_dur(rel_inds.hpc,1)),...
    'VariableNames', {'rlogP', 'rcwP', 'dprim', 'rconf', ... 
    'numNeurons', 'numTrials', 'Age', 'TaskVersion', 'ImageDuration'});

% run behav pca
behav_metrics = table2array(rank_table.hpc(:,1:4));
[behav_coeff,behav_score, behav_latent,~,behav_exp_var,~] =...
    pca(zscore(behav_metrics), 'Centered', true, 'VariableWeights', 'Variance');
loadings = (corr(behav_score, behav_metrics(:,1:size(behav_score,2))))';
exp_var_1st_pc = behav_exp_var(1);
exp_normvar_1st_pc = behav_exp_var(1)/(100/size(behav_exp_var,1));
% save in rank table
rank_table.hpc.behav_pc = behav_score(:,1);
rank_table.hpc.rank_behav_pc = tiedrank(behav_score(:,1));

% amygdala
% convert behav to rank and compute pc 
rank_table.amy = table( ... 
    tiedrank(behav_dat.session_1.rlogP(rel_inds.amy)),...
    tiedrank(behav_dat.session_1.rcwP(rel_inds.amy)),...
    tiedrank(behav_dat.session_1.dprim(rel_inds.amy)),...
    tiedrank(behav_dat.session_1.rConf(rel_inds.amy)),...
    tiedrank(neuron_nums(rel_inds.amy,1)),...
    tiedrank(trial_nums(rel_inds.amy,1)),...
    (behav_dat.session_1.Age(rel_inds.amy)),...
    categorical(spike_summary.task_vars(rel_inds.amy,1)),...
    categorical(spike_summary.image_dur(rel_inds.amy,1)),...
    'VariableNames', {'rlogP', 'rcwP', 'dprim', 'rconf', ... 
    'numNeurons', 'numTrials', 'Age', 'TaskVersion', 'ImageDuration'});

% run behav pca
behav_metrics = table2array(rank_table.amy(:,1:4));
[behav_coeff,behav_score, behav_latent,~,behav_exp_var,~] =...
    pca(zscore(behav_metrics), 'Centered', true, 'VariableWeights', 'Variance');
loadings = (corr(behav_score, behav_metrics(:,1:size(behav_score,2))))';
exp_var_1st_pc = behav_exp_var(1);
exp_normvar_1st_pc = behav_exp_var(1)/(100/size(behav_exp_var,1));
% save in rank table
rank_table.amy.behav_pc = behav_score(:,1);
rank_table.amy.rank_behav_pc = tiedrank(behav_score(:,1));

%% Fig 2, early vs late layer correlations
rel_ind = pls_results.pe.hippocampus.late_layers > 0;

scatter(pls_results.pe.hippocampus.late_layers(rel_ind), ... 
        pls_results.pe.hippocampus.early_layers(rel_ind)); 
early = pls_results.pe.hippocampus.early_layers(rel_ind); 
late = pls_results.pe.hippocampus.late_layers(rel_ind);
or_corr = corrcoef(tiedrank(late), tiedrank(early)); 
or_corr = or_corr(1,2);

n = size(early,1);
rep = 10000;
for i = 1:rep
    r = randperm(n);
    c = corrcoef(tiedrank(late(r)), tiedrank(early));
    perm_corr(i) = c(1,2);
end

perm_pval = sum(perm_corr>or_corr)/rep;
xlim([0,.8]);
ylim([0,.8]);
text(.4, .4, sprintf("p-val = %.04f", perm_pval));

%% Fig 3, original estimates
if compute_estimates==1
    % rank 
    % hippocampus
    late_lvlv = pls_results.pe.hippocampus.late_layers;                         %late
    rank_table.hpc.late_layers = tiedrank(late_lvlv(rel_inds.hpc));
    early_lvlv = pls_results.pe.hippocampus.early_layers;                       %early
    rank_table.hpc.early_layers = tiedrank(early_lvlv(rel_inds.hpc));
    all_lvlv = pls_results.pe.hippocampus.all_layers;                           %all
    rank_table.hpc.all_layers = tiedrank(all_lvlv(rel_inds.hpc));
    rate_lvlv = pls_results.spike_rate.hpc.late_layers;                         %rate
    rank_table.hpc.rate_lvlv = tiedrank(rate_lvlv(rel_inds.hpc)); 
    % amygdala
    late_lvlv = pls_results.pe.hippocampus.late_layers;                         %late
    rank_table.amy.hpc_late = tiedrank(late_lvlv(rel_inds.amy));
    late_lvlv = pls_results.pe.amygdala.late_layers;                            %amy late
    rank_table.amy.amy_late = tiedrank(late_lvlv(rel_inds.amy));
    early_lvlv = pls_results.pe.amygdala.early_layers;                          %amy early
    rank_table.amy.amy_early = tiedrank(early_lvlv(rel_inds.amy));

    % get estims
    % zero order
    zero_order_model = fitlm(rank_table.hpc, ['rank_behav_pc ~ late_layers']);  %late
    orig_estims.late_layers = zero_order_model.Coefficients.Estimate(2);
    zero_order_model = fitlm(rank_table.hpc, ['rank_behav_pc ~ early_layers']); %early
    orig_estims.early_layers = zero_order_model.Coefficients.Estimate(2);
    zero_order_model = fitlm(rank_table.hpc, ['rank_behav_pc ~ all_layers']);   %all
    orig_estims.all_layers = zero_order_model.Coefficients.Estimate(2);

    % controls
    % confounders
    perf_cont_control = fitlm(rank_table.hpc, ['rank_behav_pc ~  ' ...
        'numNeurons + numTrials + Age + TaskVersion + ImageDuration']);
    late_cont_control = fitlm(rank_table.hpc, ['late_layers ~  ' ...
        'numNeurons + numTrials + Age + TaskVersion + ImageDuration']);
    orig_estims.controls.confounder = (corr(tiedrank(perf_cont_control.Residuals.Raw),...
        tiedrank(late_cont_control.Residuals.Raw)));
    % against early 
    perf_cont_early = fitlm(rank_table.hpc, ['rank_behav_pc ~  early_layers']);
    late_cont_early = fitlm(rank_table.hpc, ['late_layers ~  early_layers']);
    orig_estims.controls.early = (corr(tiedrank(perf_cont_early.Residuals.Raw),...
        tiedrank(late_cont_early.Residuals.Raw)));
    % against all
    perf_cont_all = fitlm(rank_table.hpc, ['rank_behav_pc ~  all_layers']);
    late_cont_all = fitlm(rank_table.hpc, ['late_layers ~  all_layers']);
    orig_estims.controls.all = (corr(tiedrank(perf_cont_all.Residuals.Raw),...
        tiedrank(late_cont_all.Residuals.Raw)));
    % rate 
    perf_cont_rate = fitlm(rank_table.hpc, ['rank_behav_pc ~  rate_lvlv']);
    late_cont_rate = fitlm(rank_table.hpc, ['late_layers ~  rate_lvlv']);
    orig_estims.controls.rate = (corr(tiedrank(perf_cont_rate.Residuals.Raw),...
        tiedrank(late_cont_rate.Residuals.Raw)));
    % amygdala
    perf_cont_amy = fitlm(rank_table.amy, ['rank_behav_pc ~ amy_late + amy_early']);
    pe_cont_amy = fitlm(rank_table.amy, ['hpc_late ~ amy_late + amy_early']);
    orig_estims.controls.amygdala = (corr(tiedrank(perf_cont_amy.Residuals.Raw),...
        tiedrank(pe_cont_amy.Residuals.Raw)));
    % amygdala models
    amy_model = fitlm(rank_table.amy, ['rank_behav_pc ~ amy_early']);
    orig_estims.amygdala.early_layers = amy_model.Coefficients.Estimate(2);
    amy_model = fitlm(rank_table.amy, ['rank_behav_pc ~ amy_late']);
    orig_estims.amygdala.late_layers = amy_model.Coefficients.Estimate(2);

    % reverse models
    % early 
    perf_cont_late = fitlm(rank_table.hpc, ['rank_behav_pc ~  late_layers']);
    early_cont_late = fitlm(rank_table.hpc, ['early_layers ~  late_layers']);
    orig_estims.reverse.early = (corr(tiedrank(perf_cont_late.Residuals.Raw),...
        tiedrank(early_cont_late.Residuals.Raw)));
    % all 
    perf_cont_late = fitlm(rank_table.hpc, ['rank_behav_pc ~  late_layers']);
    all_cont_late = fitlm(rank_table.hpc, ['all_layers ~  late_layers']);
    orig_estims.reverse.all = (corr(tiedrank(perf_cont_late.Residuals.Raw),...
        tiedrank(all_cont_late.Residuals.Raw)));

    save([perm_dir + 'orig_estims.mat'], 'orig_estims')
end

%% Fig 3 compute stats for permutation runs
if compute_permstats == 1 
    n_perm = size(permutations.pe.hippocampus.all_layers,2);
    for p = 1:n_perm
        % rank 
        % hpc
        late_lvlv = permutations.pe.hippocampus.late_layers(:,p);               %late
        rank_table.hpc.late_layers = tiedrank(late_lvlv(rel_inds.hpc));
        early_lvlv = permutations.pe.hippocampus.early_layers(:,p);             %early
        rank_table.hpc.early_layers = tiedrank(early_lvlv(rel_inds.hpc));
        all_lvlv = permutations.pe.hippocampus.all_layers(:,p);                 %all
        rank_table.hpc.all_layers = tiedrank(all_lvlv(rel_inds.hpc));
        rate_lvlv = permutations.spike_rate.hpc.late_layers(:,p);               %rate
        rank_table.hpc.rate_lvlv = tiedrank(rate_lvlv(rel_inds.hpc)); 
        % amygdala
        late_lvlv = permutations.pe.hippocampus.late_layers(:,p);               %late
        rank_table.amy.hpc_late = tiedrank(late_lvlv(rel_inds.amy));
        late_lvlv = permutations.pe.amygdala.late_layers(:,p);                  %amy late
        rank_table.amy.amy_late = tiedrank(late_lvlv(rel_inds.amy));
        early_lvlv = permutations.pe.amygdala.early_layers(:,p);                %amy early
        rank_table.amy.amy_early = tiedrank(early_lvlv(rel_inds.amy));

        % permutation estimates
        % zero-order
        zero_order_model = fitlm(rank_table.hpc, ['rank_behav_pc ~ late_layers']);
        perm_estims.late_layers(p) = zero_order_model.Coefficients.Estimate(2);
        zero_order_model = fitlm(rank_table.hpc, ['rank_behav_pc ~ early_layers']);
        perm_estims.early_layers(p) = zero_order_model.Coefficients.Estimate(2);
        zero_order_model = fitlm(rank_table.hpc, ['rank_behav_pc ~ all_layers']);
        perm_estims.all_layers(p) = zero_order_model.Coefficients.Estimate(2);

        % controls
        % late against confounders
        perf_cont_control = fitlm(rank_table.hpc, ['rank_behav_pc ~  ' ...
            'numNeurons + numTrials + Age + TaskVersion + ImageDuration']);
        late_cont_control = fitlm(rank_table.hpc, ['late_layers ~  ' ...
            'numNeurons + numTrials + Age + TaskVersion + ImageDuration']);
        perm_estims.controls.confounder(p) = (corr(tiedrank(perf_cont_control.Residuals.Raw),...
            tiedrank(late_cont_control.Residuals.Raw)));
        % against early 
        perf_cont_early = fitlm(rank_table.hpc, ['rank_behav_pc ~  early_layers']);
        late_cont_early = fitlm(rank_table.hpc, ['late_layers ~  early_layers']);
        perm_estims.controls.early(p) = (corr(tiedrank(perf_cont_early.Residuals.Raw),...
            tiedrank(late_cont_early.Residuals.Raw)));
        % against all layers 
        perf_cont_all = fitlm(rank_table.hpc, ['rank_behav_pc ~  all_layers']);
        late_cont_all = fitlm(rank_table.hpc, ['late_layers ~  all_layers']);
        perm_estims.controls.all(p) = (corr(tiedrank(perf_cont_all.Residuals.Raw),...
            tiedrank(late_cont_all.Residuals.Raw)));
        % rate 
        perf_cont_rate = fitlm(rank_table.hpc, ['rank_behav_pc ~  rate_lvlv']);
        late_cont_rate = fitlm(rank_table.hpc, ['late_layers ~  rate_lvlv']);
        perm_estims.controls.rate(p) = (corr(tiedrank(perf_cont_rate.Residuals.Raw),...
            tiedrank(late_cont_rate.Residuals.Raw)));
        % amygdala
        perf_cont_amy = fitlm(rank_table.amy, ['rank_behav_pc ~ amy_late + amy_early']);
        pe_cont_amy = fitlm(rank_table.amy, ['hpc_late ~ amy_late + amy_early']);
        perm_estims.controls.amygdala(p) = (corr(tiedrank(perf_cont_amy.Residuals.Raw),...
            tiedrank(pe_cont_amy.Residuals.Raw)));

        %amygdala models
        amy_model = fitlm(rank_table.amy, ['rank_behav_pc ~ amy_early']);
        perm_estims.amygdala.early_layers(p) = amy_model.Coefficients.Estimate(2);
        amy_model = fitlm(rank_table.amy, ['rank_behav_pc ~ amy_late']);
        perm_estims.amygdala.late_layers(p) = amy_model.Coefficients.Estimate(2); 

        % reverse controls
        % early 
        perf_cont_late = fitlm(rank_table.hpc, ['rank_behav_pc ~  late_layers']);
        early_cont_late = fitlm(rank_table.hpc, ['early_layers ~  late_layers']);
        perm_estims.reverse.early(p) = (corr(tiedrank(perf_cont_late.Residuals.Raw),...
            tiedrank(early_cont_late.Residuals.Raw)));
        % all 
        perf_cont_late = fitlm(rank_table.hpc, ['rank_behav_pc ~  late_layers']);
        all_cont_late = fitlm(rank_table.hpc, ['all_layers ~  late_layers']);
        perm_estims.reverse.all(p) = (corr(tiedrank(perf_cont_late.Residuals.Raw),...
            tiedrank(all_cont_late.Residuals.Raw)));
    end
    save([perm_dir + 'permutation_estims.mat'], 'perm_estims')
end

%% plot permutation results
% load perm results 

load([perm_dir + 'permutation_estims_or.mat']);
load([perm_dir + 'orig_estims_or.mat']);

% estim p-values and plot null distributions
nbins = 15;
figure

% early layers
a = subplot(3,4,1); 
perm_pval = sum(perm_estims.early_layers>.32) ... 
    /length(perm_estims.early_layers);
perm_pvals.early_layers = perm_pval; 
histogram(perm_estims.early_layers, nbins, 'FaceAlpha', 0.5);
xlabel("beta estimates"); 
title("Early layers"); 
xline(orig_estims.early_layers, "Linewidth", 2)
xxs = xlim(); yys = ylim();
text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2), sprintf("p-val = %.03f", round(perm_pval*1000)/1000));

% late layers
b = subplot(3,4,2);
perm_pval = sum(perm_estims.late_layers>.54) ...  
    /length(perm_estims.late_layers);
perm_pvals.late_layers = perm_pval;
histogram(perm_estims.late_layers, nbins, 'FaceAlpha', 0.5);
xlabel("beta estimates"); 
title("Late layers"); 
xline(orig_estims.late_layers, "Linewidth", 2)
xxs = xlim(); yys = ylim();
text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2), sprintf("p-val = %.03f", round(perm_pval*1000)/1000));

% all layers
c = subplot(3,4,3);
perm_pval = sum(perm_estims.all_layers>.42) ...  
    /length(perm_estims.all_layers);
perm_pvals.all_layers = perm_pval;
histogram(perm_estims.all_layers, nbins, 'FaceAlpha', 0.5);
xlabel("beta estimates"); 
title("All layers"); 
xline(orig_estims.all_layers, "Linewidth", 2)
xxs = xlim(); yys = ylim();
text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2), sprintf("p-val = %.03f", round(perm_pval*1000)/1000));

% control confounders
d = subplot(3,4,4);
perm_pval = sum(perm_estims.controls.confounder>.64) ... 
    /length(perm_estims.controls.confounder);
perm_pvals.controls.confounder = perm_pval;
histogram(perm_estims.controls.confounder, nbins, 'FaceAlpha', 0.5);
xlabel("partial eta"); 
title("Control Confounders","Num trial & neuron, age, task, dur"); 
xline(orig_estims.controls.confounder, "Linewidth", 2)
xxs = xlim(); yys = ylim();
text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2), sprintf("p-val = %.03f", round(perm_pval*1000)/1000));
display(max(perm_estims.controls.confounder));

% control early
e = subplot(3,4,5);
perm_pval = sum(perm_estims.controls.early>.45) ... 
    /length(perm_estims.controls.early);
perm_pvals.controls.early = perm_pval;
histogram(perm_estims.controls.early, nbins, 'FaceAlpha', 0.5);
xlabel("partial eta"); 
title("Control early"); 
xline(orig_estims.controls.early, "Linewidth", 2)
xxs = xlim(); yys = ylim();
text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2), sprintf("p-val = %.03f", round(perm_pval*1000)/1000));

% control rate
f = subplot(3,4,6);
perm_pval = sum(perm_estims.controls.rate>.39) ... 
    /length(perm_estims.controls.rate);
perm_pvals.controls.rate = perm_pval;
histogram(perm_estims.controls.rate, nbins, 'FaceAlpha', 0.5);
xlabel("partial eta"); 
title("Control spike rate"); 
xline(orig_estims.controls.rate, "Linewidth", 2)
xxs = xlim(); yys = ylim();
text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2), sprintf("p-val = %.03f", round(perm_pval*1000)/1000));

% amygdala early layers
g = subplot(3,4,7);
perm_pval = sum(perm_estims.amygdala.early_layers>-.16)...
    /length(perm_estims.amygdala.early_layers);
perm_pvals.amygdala.early_layers = perm_pval;
histogram(perm_estims.amygdala.early_layers, nbins, 'FaceAlpha', 0.5);
xlabel("partial eta"); 
title("Amygdala & early layers"); 
xline(orig_estims.amygdala.early_layers, "Linewidth", 2)
xxs = xlim(); yys = ylim();
text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2), sprintf("p-val = %.03f", round(perm_pval*1000)/1000));

% amygdala late layers
h = subplot(3,4,8);
perm_pval = sum(perm_estims.amygdala.late_layers>-.17) ... 
    /length(perm_estims.amygdala.late_layers);
perm_pvals.amygdala.late_layers = perm_pval;
histogram(perm_estims.amygdala.late_layers, nbins, 'FaceAlpha', 0.5);
xlabel("partial eta"); 
title("Amygdala & late layers"); 
xline(orig_estims.amygdala.early_layers, "Linewidth", 2)
xxs = xlim(); yys = ylim();
text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2), sprintf("p-val = %.03f", round(perm_pval*1000)/1000));

% control amygdala
i = subplot(3,4,9);
perm_pval = sum(perm_estims.controls.amygdala>.57) ... 
    /length(perm_estims.controls.amygdala);
perm_pvals.controls.amygdala = perm_pval;
histogram(perm_estims.controls.amygdala, nbins, 'FaceAlpha', 0.5);
xlabel("partial eta"); 
title("HPC controlled for Amygdala"); 
xline(orig_estims.controls.amygdala, "Linewidth", 2)
xxs = xlim(); yys = ylim();
text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2), sprintf("p-val = %.03f", round(perm_pval*1000)/1000));

% additional controls which are not in the main plot 
% control for all layers
j = subplot(3,4,10);
perm_pval = sum(perm_estims.controls.all>.38) ...  
    /length(perm_estims.controls.all);
perm_pvals.controls.all = perm_pval;
histogram(perm_estims.controls.all, nbins, 'FaceAlpha', 0.5);
xlabel("partial eta"); 
title("Late controlled for All"); 
xline(orig_estims.controls.all, "Linewidth", 2)
xxs = xlim(); yys = ylim();
text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2), sprintf("p-val = %.04f", round(perm_pval*10000)/10000));

% control early for late 
k = subplot(3,4,11);
perm_pval = sum(perm_estims.reverse.early>-.02) ... 
    /length(perm_estims.reverse.early);
perm_pvals.reverse.early = perm_pval;
histogram(perm_estims.reverse.early, nbins, 'FaceAlpha', 0.5);
xlabel("partial eta"); 
title("R: Early controlled for Late"); 
xline(orig_estims.reverse.early, "Linewidth", 2)
xxs = xlim(); yys = ylim();
text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2), sprintf("p-val = %.03f", round(perm_pval*1000)/1000));

% control all for late
l = subplot(3,4,12);
perm_pval = sum(perm_estims.reverse.all>-.11) ... 
    /length(perm_estims.reverse.all);
perm_pvals.reverse.all = perm_pval;
histogram(perm_estims.reverse.all, nbins, 'FaceAlpha', 0.5);
xlabel("partial eta"); 
title("R: All controlled for Late"); 
xline(orig_estims.reverse.all, "Linewidth", 2)
xxs = xlim(); yys = ylim();
text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2), sprintf("p-val = %.03f", round(perm_pval*1000)/1000));


linkaxes([a, b, c, d, e, f, g, h, i, j, k, l], 'x');



