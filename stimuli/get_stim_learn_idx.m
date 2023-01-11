function stim_learn_idx = get_stim_learn_idx(all_HMAX_info, fulltable)
    %% get indices of learned stimuli
    %% ordered as they were presented to subjects
    for il = 1:size(fulltable,1)
        stim_learn_idx(il) = find(all_HMAX_info.StimulusCode==fulltable.stimuliLearn(il));
    end
end