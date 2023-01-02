function ori_estim = get_ori_wise_estim(c1, c1Scale, orientations)
%ORI_ESTIM Get orientation-wise entropy and median 
    for k=1:length(c1Scale)-1
        for j=1:length(orientations)
                ori_estim.c1_median(k, j)=median(nonzeros(c1{1}{1, k}(:, :, j)));
                ori_estim.c1_sd(k, j)=std(nonzeros(c1{1}{1, k}(:, :, j)));
                ori_estim.c1_median_entropy(k, j) = median(nonzeros(entropyfilt((c1{1}{1, k}(:, :, j)))));
                ori_estim.c1_entropy(k,j) = entropy(c1{1}{1, k}(:,:,j));
        end
    end
end

