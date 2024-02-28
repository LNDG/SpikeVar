function avg_behav_dat = a_extract_behav_info(fulltable_learn, fulltable_recog, pat_ids, behav_dir)
%%extract_behav_info Takes as input the fulltable_learn and fulltable_recog
%   and extracts all the relevant behavioral information for the behavioral
%   models.

for pat_id = 1:length(pat_ids)
    c_id = pat_ids(pat_id);
    % get relevant indices
    % set of neurons (all, hippocampus <3  or amygdala >2)
    for i_ses = 1:2
        % recog
        all_ids = find(fulltable_recog.Participant==c_id &...
            fulltable_recog.Session==i_ses);
        if ~isempty(all_ids)
            rel_ids  = all_ids(1:fulltable_recog.NumTrials(all_ids(1)));
            rel_RT = (fulltable_recog.RecogRTAll(rel_ids))/1000;
            rel_RS = ((fulltable_recog.RecogRTAll(rel_ids))/1000).^-1;
            med_all_Recog_RT(pat_id,i_ses) = median(rel_RT);
            med_all_Recog_RS(pat_id,i_ses) = median(rel_RS);
            % get accuracy
            resp_given = fulltable_recog.RecogResponseAll(rel_ids);
            resp_given(resp_given<=3) = 0; resp_given(resp_given>3) = 1;
            stim_cat = fulltable_recog.RecogOldNew(rel_ids);
            resp_eval = zeros(size(stim_cat)); resp_eval(resp_given==stim_cat) = 1;
            recog_all_accuracy(pat_id,i_ses) = nanmean(resp_eval);
            recog_all_logit_accuracy(pat_id,i_ses) = log(recog_all_accuracy(pat_id,i_ses)/...
                (1-recog_all_accuracy(pat_id,i_ses)));
            % confidence
            recog_all_conf = fulltable_recog.RecogResponseAll(rel_ids);
            recog_all_conf(recog_all_conf<=3) = abs(recog_all_conf(recog_all_conf<=3)-4);
            recog_all_conf(recog_all_conf>3) = recog_all_conf(recog_all_conf>3)-3;
            recog_avg_conf(pat_id,i_ses) = mean(recog_all_conf);
            % confidence separately for correct and incorrect
            recog_avg_corr_conf(pat_id,i_ses) = mean(recog_all_conf(resp_eval==1));
            recog_avg_incorr_conf(pat_id,i_ses) = mean(recog_all_conf(resp_eval==0));
            
            % confidence weighted accuracy
            recog_conf_wght_acc(pat_id, i_ses) = nanmean(recog_all_conf.*resp_eval);
            
            
            % check within person correlation between confidence and RT and
            % RS
            [conf_rt_corr(pat_id,i_ses),p] = corr((fulltable_recog.RecogRTAll(rel_ids))/1000,...
                recog_all_conf);
            [conf_rs_corr(pat_id,i_ses),p] = corr(((fulltable_recog.RecogRTAll(rel_ids))/1000).^-1,...
                recog_all_conf);
            % to make this even easier, ectract median RT and RS for
            % confidence bins (1-3)
            for confval = 1:3
                conf_dep_med_RT(pat_id, i_ses, confval) = nanmedian(rel_RT(recog_all_conf==confval));
                conf_dep_med_RS(pat_id, i_ses, confval) = nanmedian(rel_RS(recog_all_conf==confval));
            end           
            % loop across trials and determine performance
            FAs = NaN(1,length(rel_ids)); CRs = NaN(1,length(rel_ids));
            Hs = NaN(1,length(rel_ids)); Ms = NaN(1,length(rel_ids));
            for itr = 1:length(rel_ids)
                % new = 0, old = 1 in RecogOldNew
                % HIT = old, old
                if fulltable_recog.RecogResponseAll(rel_ids(itr))>3&&...
                        fulltable_recog.RecogOldNew(rel_ids(itr))==1
                    Hs(itr) = 1;
                    % False alaram = old response, new stim
                elseif fulltable_recog.RecogResponseAll(rel_ids(itr))>3&&...
                        fulltable_recog.RecogOldNew(rel_ids(itr))==0
                    FAs(itr) = 1;
                    % Correct rejection = new, new
                elseif fulltable_recog.RecogResponseAll(rel_ids(itr))<=3&&...
                        fulltable_recog.RecogOldNew(rel_ids(itr))==0
                    CRs(itr) = 1;
                    % Miss = old stim, new response
                elseif fulltable_recog.RecogResponseAll(rel_ids(itr))<=3&&...
                        fulltable_recog.RecogOldNew(rel_ids(itr))==1
                    Ms(itr) = 1;
                end
            end
            
            Hrate(pat_id,i_ses) = nansum(Hs)/(nansum(Hs)+nansum(Ms));
            FArate(pat_id,i_ses) = nansum(FAs)/(nansum(FAs)+nansum(CRs));
            % d prime
            Dprime(pat_id,i_ses) = norminv(Hrate(pat_id,i_ses))-...
                norminv(FArate(pat_id,i_ses));
            % criterion
            Criterion(pat_id, i_ses) = -0.5*(norminv(Hrate(pat_id,i_ses))+...
                norminv(FArate(pat_id,i_ses)));
            
            age(pat_id, i_ses) = unique(fulltable_recog.Age(rel_ids));
            id_inf(pat_id, i_ses) = unique(fulltable_recog.Participant(rel_ids));
            
        else
            Hrate(pat_id,i_ses) = NaN;
            FArate(pat_id,i_ses) = NaN;
            Dprime(pat_id,i_ses) = NaN;
            Criterion(pat_id, i_ses) = NaN;
            
            recog_all_accuracy(pat_id,i_ses) = NaN;
            recog_all_logit_accuracy(pat_id,i_ses) = NaN;
            
            recog_avg_conf(pat_id,i_ses) = NaN;
            recog_avg_corr_conf(pat_id,i_ses) = NaN;
            recog_avg_incorr_conf(pat_id,i_ses) = NaN;
            recog_conf_wght_acc(pat_id, i_ses) = NaN;
            
            med_all_Recog_RT(pat_id,i_ses) = NaN;
            med_all_Recog_RS(pat_id,i_ses) = NaN;
            
            age(pat_id, i_ses) = NaN;
            id_inf(pat_id, i_ses) = NaN;
        end
        
        all_l_ids = find(fulltable_learn.Participant==c_id &...
            fulltable_learn.Session==i_ses);
        if ~isempty(all_l_ids)
            rel_ids  = all_l_ids(1:fulltable_learn.NumTrials(all_l_ids(1)));
            rel_RT = (fulltable_learn.learnRT(rel_ids))/1000;
            rel_RS = ((fulltable_learn.learnRT(rel_ids))/1000).^-1;
            med_all_Learn_RT(pat_id,i_ses) = median(rel_RT);
            med_all_Learn_RS(pat_id,i_ses) = median(rel_RS);
        else
            med_all_Learn_RT(pat_id,i_ses) = NaN;
            med_all_Learn_RS(pat_id,i_ses) = NaN;          
        end
        
        % learning
        all_ids = find(fulltable_learn.Participant==c_id &...
            fulltable_learn.Session==i_ses);
        if ~isempty(all_ids)
            rel_ids  = all_ids(1:fulltable_recog.NumTrials(all_ids(1)));
            % learn accuracy 
            learn_accuracy(pat_id, i_ses) = nanmean(fulltable_learn.learnCorrect(rel_ids)); 
        end         
    end
end
%% convert to table
avg_behav_dat = struct('session_1', array2table(NaN(42,10), 'VariableNames',...
     {'ID','Age','lRT', 'rRT', ...
     'lRS', 'rRS', 'rPerf', 'rlogP', 'rConf', 'rcwP'}), ... 
    'session_2', array2table(NaN(42,10), 'VariableNames',...
     {'ID','Age','lRT', 'rRT', ...
     'lRS', 'rRS', 'rPerf', 'rlogP', 'rConf', 'rcwP'})); 
 
for i_ses = 1:2
    field_name = sprintf('session_%d', i_ses);
    avg_behav_dat.(field_name).rRT = med_all_Recog_RT(:,i_ses);
    avg_behav_dat.(field_name).rRS = med_all_Recog_RS(:,i_ses);
    
    avg_behav_dat.(field_name).lRT = med_all_Learn_RT(:,i_ses);
    avg_behav_dat.(field_name).lRS = med_all_Learn_RS(:,i_ses);
    
    avg_behav_dat.(field_name).rPerf = recog_all_accuracy(:,i_ses);
    avg_behav_dat.(field_name).rlogP = recog_all_logit_accuracy(:,i_ses);
    
    avg_behav_dat.(field_name).rConf = recog_avg_conf(:,i_ses);
    avg_behav_dat.(field_name).rcwP = recog_conf_wght_acc(:,i_ses);
    
    % add criterion and dprime
    avg_behav_dat.(field_name).crit = Criterion(:,i_ses);
    avg_behav_dat.(field_name).dprim = Dprime(:,i_ses);
    
    % add age and ID
    avg_behav_dat.(field_name).ID = id_inf(:,i_ses);
    avg_behav_dat.(field_name).Age = age(:,i_ses);
    
    % learn accuracy 
    avg_behav_dat.(field_name).learn_accuracy = learn_accuracy(:,i_ses);
end
%% save 
save([behav_dir + 'SpikeVar_behav_info.mat'], 'avg_behav_dat');
end