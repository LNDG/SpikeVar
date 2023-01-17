function patch_estim = get_patch_wise_estim(c1,c2)
% calculate responses for each patch
    c1_mean = cellfun(@(x) mean(x,3), c1{1}, 'uni', 0); % take ori 4
    
    patch_estim.c1_median = cell2mat(cellfun(@(x) median(x(:)), c1_mean, 'uni', 0)); % take median
    patch_estim.c2_median = cellfun(@median, c2, 'uni', 0); % take median
    
    patch_estim.c1_sum = cell2mat(cellfun(@(x) sum(x(:)), c1_mean, 'uni', 0)); % take sum
    patch_estim.c2_sum  = cell2mat(cellfun(@sum, c2, 'uni', 0)); % take sum

    patch_estim.c1_sd = cell2mat(cellfun(@(x) std(x(:)), c1_mean, 'uni', 0)); % take std
    patch_estim.c2_sd = cellfun(@std, c2, 'uni', 0); % take std    
end

