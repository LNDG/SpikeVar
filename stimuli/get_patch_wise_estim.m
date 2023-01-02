function patch_estim = get_patch_wise_estim(c1,c2)
% calculate responses for each patch
    c1_ori4 = cellfun(@(x) x(:,:,4), c1{1}, 'uni', 0); % take ori 4
    
    patch_estim.c1_median = cell2mat(cellfun(@(x) median(x(:)), c1_ori4, 'uni', 0)); % take median
    patch_estim.c2_median = cellfun(@median, c2, 'uni', 0); % take median

    patch_estim.c1_sd = cell2mat(cellfun(@(x) std(x(:)), c1_ori4, 'uni', 0)); % take std
    patch_estim.c2_sd = cellfun(@std, c2, 'uni', 0); % take std
end

