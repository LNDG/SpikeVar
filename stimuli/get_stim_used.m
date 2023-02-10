function stim_used = get_stim_used(info, fulltable_learn, fulltable_recog)
%% get all stimulus that were presented during encoding and recognition
    n_images = size(info,1);
    for i = 1:n_images
       stim_learn(i) = ismember(info.stimulus_code(i), fulltable_learn.stimuliLearn); 
       stim_recog(i) = ismember(info.stimulus_code(i), fulltable_recog.stimuliRecog);
    end
    stim_used  = stim_learn + stim_recog;
end