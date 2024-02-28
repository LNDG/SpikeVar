function d_lm_layer_vs_perf(rank_table)
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

    % plot
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
        scatter(table2array(rank_table.hpc(:,imod)), rank_table.hpc.rank_behav_pc,140,...
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
end

