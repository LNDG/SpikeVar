function all_stim = get_stim_used(all_HMAX_info, fulltable_learn, fulltable_recog)
%% get stimulus that were presented during encoding
    n_images = size(all_HMAX_info,1) 
    for i = 1:n_images
       stim_learn(i) = ismember(all_HMAX_info.StimulusCode(i), fulltable_learn.stimuliLearn); % refactor
       stim_recog(i) = ismember(all_HMAX_info.StimulusCode(i), fulltable_recog.stimuliRecog);
    end
    all_stim  = stim_learn + stim_recog;
end