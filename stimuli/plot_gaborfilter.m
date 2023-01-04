function plot_gaborfilter(filter, filtersize)
% function to plot the gaborfilter
    filter = filter(filter~=0);
    filter = reshape(filter, [filtersize, filtersize]);
    imagesc(filter);
    colorbar();
    colormap(gray);
end