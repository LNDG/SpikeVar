function plot_lvlv_corrs(lvlv_corrs, lvlv_corr_lims, plot_ids, plot_title)
%% Plot latent correlation including CIs
    figure
    scatter(1:length(plot_ids), lvlv_corrs(plot_ids), 100,...
        'filled', 'MarkerFaceAlpha', .8, 'MarkerFaceColor', [.5, .3, .7]);
    hold on
    for i_pat = 1:length(plot_ids)
       plot([(i_pat) (i_pat)], ...
           lvlv_corr_lims(plot_ids(i_pat),:), 'k') 
    end
    ylim([0 1])
    yticks([0:.2:1])
    xs = xlim;
    plot(xs, [0 0], '--', 'Color', [.7 .7 .7])
    xlabel('Participants')
    ylabel('Absolute latent correlation')
    title(plot_title)
end 