function plot_raincloud_stim_weights(stim_weights, layers, plot_ids)
%% Raincloud plots of absolute stimulus weights
    xdat = [];
    ydat = [];
    % First for Hippocampus PE and late layers
    % re-arrange data first
    for sub_id = 1:size(stim_weights,1)
        if ~isempty(stim_weights{sub_id,1})
            all_abs_stim_weights(:,sub_id) = (abs(stim_weights{sub_id,1}(:,1)));
        else
            all_abs_stim_weights(:,sub_id) = NaN(size(stim_weights{1,1},1),1);
        end
    end

    for lay = 1:size(stim_weights{1},1)
        abs_stim_wghts{lay,1} = (all_abs_stim_weights(lay,plot_ids))';
    end

    cl = [.3 .4 .7];
    figure
    hold on
    h   = rm_raincloud(abs_stim_wghts, cl);
    set(gca, 'XLim', [-0.4 1]);
    title(['Absolute stimulus weights' newline 'across ' layers ' layers']);
    
    if layers == "all"
        n_cols = [ 203, 163, 40; 35, 12, 15;  203, 163, 40; 35, 12, 15;...
            242, 220, 93;...
            242, 163, 89;...
            219, 144, 101;...
            192, 74, 66;...
            164, 3, 31;...
            242, 220, 93;...
            242, 163, 89;...
            219, 144, 101;...
            192, 74, 66;...
            164, 3, 31;...
             242, 220, 93;...
            242, 163, 89;...
            219, 144, 101;...
            192, 74, 66;...
            164, 3, 31]/256;
    else
        n_cols = [135, 35, 49;...
            135, 35, 49;...
            [255, 200, 87;...
            233, 114, 76;...
            197, 40, 61];...
            [255, 200, 87;...
            233, 114, 76;...
            197, 40, 61];...
            [255, 200, 87;...
            233, 114, 76;...
            197, 40, 61]]/255;        
    end

    % change colors
    for lay = 1:size(stim_weights{1},1)
    h.p{lay, 1}.FaceColor         = n_cols(lay,:);
    h.s{lay, 1}.MarkerFaceColor   = n_cols(lay,:);
    h.m(lay, 1).MarkerEdgeColor   = 'none';
    h.m(lay, 1).MarkerFaceColor   = n_cols(lay,:);
        if lay<size(stim_weights{1},1)
            h.l(lay, 1).Color              = [.7 .7 .7];
        end
        
        if layers == "all"
            if lay == 2 || lay == 4 || lay == 9  || lay == 14 
                h.l(lay,1).XData = [];
                h.l(lay,1).YData = [];
            end
        else
            if lay == 2 || lay == 5 || lay == 8
                h.l(lay,1).XData = [];
                h.l(lay,1).YData = [];
            end
        end

    end

    ylabel('Image features')
    
    if layers == "late"
        yticklabels({'VGG 5 Num Z', 'VGG 4 Num Z',...
        'VGG 3 Num Z','VGG 5 SD', 'VGG 4SD',...
        'VGG 3 SD', 'VGG 5 Sum', 'VGG 4 Sum',...
        'VGG 3 Sum','C2 SD', 'C2 Sum'})
    elseif layers == "early"
        yticklabels({'VGG 5 Num Z', 'VGG 4 Num Z',...
        'VGG 3 Num Z','VGG 5 SD', 'VGG 4SD',...
        'VGG 3 SD', 'VGG 5 Sum', 'VGG 4 Sum',...
        'VGG 3 Sum','C1 SD', 'C1 Sum'})
    else 
        yticklabels({ 'VGG 5 Z Num', 'VGG 4 Z Num',...
        'VGG 3 Z Num','VGG 2 Z Num','VGG 1 Z Num',...
        'VGG 5 SD', 'VGG 4SD',...
        'VGG 3 SD', 'VGG 2 SD', 'VGG 1 SD', 'VGG 5 Sum', 'VGG 4 Sum',...
        'VGG 3 Sum','VGG 2 Sum','VGG 1 Sum',...
        'C2 SD', 'C1 SD','C2 Sum', 'C1 Sum'})
    end 
    
    xlabel('Absolute stimulus weight')



end
