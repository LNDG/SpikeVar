function plot_lvlv_corrs_all_layers(pls_results)
%% Supplement: all layer Hippocampus PE
    plot_ids = find(pls_results.lvlv_corrs.all_layers{1}>0);
    figure
    subplot(1,3,1)
    scatter(1:length(plot_ids),pls_results.lvlv_corrs.early_layers{1}(plot_ids), 100,...
        'filled', 'MarkerFaceAlpha', .8, 'MarkerFaceColor', [.5, .3, .7]);
    hold on
    for i_pat = 1:length(plot_ids)
       plot([(i_pat) (i_pat)],...
           pls_results.lvlv_corr_lims.early_layers{1}(plot_ids(i_pat),:), 'k') 
    end
    ylim([-.2 1])
    yticks([0:.2:1])
    xs =xlim;
    plot(xs, [0 0], '--', 'Color', [.7 .7 .7])
    xlabel('Participants')
    ylabel('Absolute latent correlation')
    title('Hippocampus PE vs early layers')


    subplot(1,3,2)
    scatter(1:length(plot_ids),pls_results.lvlv_corrs.late_layers{1}(plot_ids), 100,...
        'filled', 'MarkerFaceAlpha', .8, 'MarkerFaceColor', [.5, .3, .7]);
    hold on
    for i_pat = 1:length(plot_ids)
       plot([(i_pat) (i_pat)],...
           pls_results.lvlv_corr_lims.late_layers{1}(plot_ids(i_pat),:), 'k') 
    end
    ylim([-.2 1])
    yticks([0:.2:1])
    xs =xlim;
    plot(xs, [0 0], '--', 'Color', [.7 .7 .7])
    xlabel('Participants')
    ylabel('Absolute latent correlation')
    title('Hippocampus PE vs late layers')


    subplot(1,3,3)
    scatter(1:length(plot_ids),pls_results.lvlv_corrs.all_layers{1}(plot_ids), 100,...
        'filled', 'MarkerFaceAlpha', .8, 'MarkerFaceColor', [.5, .3, .7]);
    hold on
    for i_pat = 1:length(plot_ids)
       plot([(i_pat) (i_pat)],...
           pls_results.lvlv_corr_lims.all_layers{1}(plot_ids(i_pat),:), 'k') 
    end
    ylim([-.2 1])
    yticks([0:.2:1])
    xs =xlim;
    plot(xs, [0 0], '--', 'Color', [.7 .7 .7])
    xlabel('Participants')
    ylabel('Absolute latent correlation')
    title('Hippocampus PE vs all layers')
end