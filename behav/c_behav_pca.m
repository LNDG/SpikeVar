function rank_table = c_behav_pca(rank_table, rel_inds)
%%b_behav_pca Performs PCA on behavioral measures 
    behav_metrics = table2array(rank_table.hpc(:,7:10));
    [behav_coeff,behav_score,...
        behav_latent,~,behav_exp_var,~] =...
        pca(zscore(behav_metrics), 'Centered', true, 'VariableWeights', 'Variance');
    loadings = (corr(behav_score,...
        behav_metrics(:,1:size(behav_score,2))))';
    exp_var_1st_pc = behav_exp_var(1);
    exp_normvar_1st_pc = behav_exp_var(1)/...
        (100/size(behav_exp_var,1));

%    spss_scores = [-.49096 -1.03220 -.95593 -1.34061 -.99076 1.58337 .08288,...
%        -.88392 -1.54987 .25994 -.72726 .07078 -1.24414 .23109 .13330 .96957,...
%        1.60498 -1.65079 -.11430 .85260 .91304 .15395 -1.07072 -.75862 1.30807,...
%        .69051 .60916 -.16001 -.79968 1.65848 .52723 .34587 1.77494];
%    corr(behav_score(:,1),spss_scores');

    rank_table.hpc.behav_pc = behav_score(:,1);
    rank_table.hpc.rank_behav_pc = tiedrank(behav_score(:,1));

    rank_table.amyg.behav_pc = behav_score(rel_inds.amyg,1);
    rank_table.amyg.rank_behav_pc = tiedrank(behav_score(rel_inds.amyg,1));
end