function [p, partial_etas, data] = e_control_analysis_thresholds(rank_table)
%e_control_analysis Runs the following control analysis: 

%% Zero order model
% Linear models: each layer against performance PC
    for imod = 1:3
        PE_perf_model_ranks{imod} = fitlm(rank_table.hpc, ['rank_behav_pc ~  ' ...
            rank_table.hpc.Properties.VariableNames{imod}]);
        PE_mod_eval.estims.behav_pc(imod) = PE_perf_model_ranks{imod}.Coefficients.Estimate(2);
        PE_mod_eval.tval.behav_pc(imod) = PE_perf_model_ranks{imod}.Coefficients.tStat(2);
        PE_mod_eval.r2.behav_pc(imod) = PE_perf_model_ranks{imod}.Rsquared.Adjusted;
        PE_mod_eval.p.behav_pc(imod) = PE_perf_model_ranks{imod}.Coefficients.pValue(2);
        temp = coefCI(PE_perf_model_ranks{imod}); PE_mod_eval.CI.behav_pc(imod,:) = temp(2,:);
    end

%% Contrast layer-specific explanatory power
% model performance as a function of all layers first
PE_perf_model_early_late_ranks = fitlm(rank_table.hpc, ['rank_behav_pc ~  ' ...
    'PE_late + PE_early ']);
PE_perf_model_all_late_ranks = fitlm(rank_table.hpc, ['rank_behav_pc ~  ' ...
    'PE_late + PE_all ']);
% Late layers survive as only predictor in both cases

% get partial etas
PE_perf_cont_early_mod_ = fitlm(rank_table.hpc, ['rank_behav_pc ~  ' ...
    'PE_early ']);
PE_perf_cont_all_mod_ = fitlm(rank_table.hpc, ['rank_behav_pc ~  ' ...
    'PE_all ']);
PE_pe_cont_early_mod_ = fitlm(rank_table.hpc, ['PE_late ~  ' ...
    'PE_early ']);
PE_pe_cont_all_mod_ = fitlm(rank_table.hpc, ['PE_late ~  ' ...
    'PE_all ']);
PE_compare.eta.early_vs_late = (corr(tiedrank(PE_perf_cont_early_mod_.Residuals.Raw),...
    tiedrank(PE_pe_cont_early_mod_.Residuals.Raw)));
PE_compare.eta.all_vs_late = (corr(tiedrank(PE_perf_cont_all_mod_.Residuals.Raw),...
    tiedrank(PE_pe_cont_all_mod_.Residuals.Raw)));

% get correct p-value
% first transform r to t
% we need the "real" n
n_actual = size(rank_table.hpc,1)-...
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
PE_perf_model_spike_late_ranks = fitlm(rank_table.hpc, ['rank_behav_pc ~  ' ...
    'PE_late + SpikeNum_late ']);
% late layer PE survives

% Get partial etas
PE_perf_cont_spiken_mod = fitlm(rank_table.hpc, ['rank_behav_pc ~  ' ...
    'SpikeNum_late']);
PE_pe_cont_spiken_mod = fitlm(rank_table.hpc, ['PE_late ~  ' ...
    'SpikeNum_late ']);

PE_compare.eta.spikenum_vs_PElate = (corr(tiedrank(PE_perf_cont_spiken_mod.Residuals.Raw),...
    tiedrank(PE_pe_cont_spiken_mod.Residuals.Raw)));

% get correct p-value
% first transform r to t
% we need the "real" n
n_actual = size(rank_table.hpc,1)-...
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

PE_perf_model_amyg_late_ranks = fitlm(rank_table.amyg, ['rank_behav_pc ~  ' ...
    'PE_late + Amyg_PE_late + Amyg_PE_early']);
% late layer PE survives

% Get partial etas
PE_perf_cont_amyg_mod = fitlm(rank_table.amyg, ['rank_behav_pc ~  ' ...
    'Amyg_PE_late + Amyg_PE_early']);
PE_pe_cont_amyg_mod = fitlm(rank_table.amyg, ['PE_late ~  ' ...
    'Amyg_PE_late + Amyg_PE_early']);

PE_compare.eta.AmygLate_vs_PElate = (corr(tiedrank(PE_perf_cont_amyg_mod.Residuals.Raw),...
    tiedrank(PE_pe_cont_amyg_mod.Residuals.Raw)));

% get correct p-value
% first transform r to t
% we need the "real" n
n_actual = size(rank_table.amyg,1)-...
    (size(PE_perf_model_amyg_late_ranks.Coefficients,1)-2);
PE_compare.t.AmygLate_vs_PElate = r_to_t(PE_compare.eta.AmygLate_vs_PElate, n_actual);
PE_compare.p.AmygLate_vs_PElate = t_to_p(PE_compare.t.AmygLate_vs_PElate, n_actual);

% Boostrap CI
mycorr = @(x1,x2) corr(x1,x2,'Type','Pearson');
nI = 10000;
[PE_compare.etaCI.AmygLate_vs_PElate , ~] = bootci(nI,{mycorr,tiedrank(PE_perf_cont_amyg_mod.Residuals.Raw),...
    tiedrank(PE_pe_cont_amyg_mod.Residuals.Raw)});

%% control late PE for additional potential confounders
PE_perf_model_late_controls_ranks = fitlm(rank_table.hpc, ['rank_behav_pc ~  ' ...
    'PE_late +  numNeurons + numTrials + Age + TaskVersion + ImageDuration']);
% PE survives

% Get partial etas
PE_perf_cont_control_mod = fitlm(rank_table.hpc, ['rank_behav_pc ~  ' ...
    'numNeurons + numTrials + Age + TaskVersion + ImageDuration']);
PE_pe_cont_control_mod_ = fitlm(rank_table.hpc, ['PE_late ~  ' ...
    'numNeurons + numTrials + Age + TaskVersion + ImageDuration']);
PE_compare.eta.PElate_controls = (corr(tiedrank(PE_perf_cont_control_mod.Residuals.Raw),...
    tiedrank(PE_pe_cont_control_mod_.Residuals.Raw)));

% get correct p-value
% first transform r to t
% we need the "real" n
n_actual = size(rank_table.hpc,1)-...
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
% zero order
figure
subplot(2,3,1)
hold on
scatter(table2array(rank_table.hpc(:,2)), rank_table.hpc.rank_behav_pc,140,...
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


% against early
subplot(2,3,2)
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

% SPIKENUM 
subplot(2,3,3)
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

% CONFOUNDERS 
subplot(2,3,4)
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

% AMYGDALA 
subplot(2,3,5)
scatter(tiedrank(PE_perf_cont_amyg_mod.Residuals.Raw),...
    tiedrank(PE_pe_cont_amyg_mod.Residuals.Raw),140,...
    'filled', 'MarkerFaceAlpha', 1, 'MarkerFaceColor', cols.pe(2,:))
lsline
xticks([])
xlabel('PE modulation [ranks]')
title('Controlled for amygdala')
ylabel('Performance [ranks]')
xxs = xlim(); yys = ylim();

text(xxs(1)+.3*xxs(2), yys(1)+.5*yys(2),['partial eta = ' num2str(...
    round(PE_compare.eta.AmygLate_vs_PElate*1000)/1000)...
    ' p = ' num2str(round(PE_compare.p.AmygLate_vs_PElate*10000)/10000)])


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

%bar plot
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

partial_etas.zero_order = round(PE_mod_eval.estims.behav_pc(2)*1000)/1000;
p.zero_order = round(PE_mod_eval.p.behav_pc(2)*1000)/1000;
partial_etas.confounders = round(PE_compare.eta.PElate_controls*1000)/1000;
p.confounders = round(PE_compare.p.PElate_controls*10000)/10000;
partial_etas.spikenum = round(PE_compare.eta.spikenum_vs_PElate*1000)/1000;
p.spikenum = round(PE_compare.p.spikenum_vs_PElate*1000)/1000;

data.zero_order = [table2array(rank_table.hpc(:,2)), rank_table.hpc.rank_behav_pc]; 
data.confounders = [tiedrank(PE_perf_cont_control_mod.Residuals.Raw), tiedrank(PE_pe_cont_control_mod_.Residuals.Raw)]; 
data.spikenum = [tiedrank(PE_perf_cont_spiken_mod.Residuals.Raw), tiedrank(PE_pe_cont_spiken_mod.Residuals.Raw)]; 

end

