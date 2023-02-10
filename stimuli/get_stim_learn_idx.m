function stim_learn_idx = get_stim_learn_idx(info, fulltable)
    %% get indices of learned stimuli
    %% ordered as they were presented to subjects
    for il = 1:size(fulltable,1)
        stim_learn_idx(il) = find(info.stimulus_code==fulltable.stimuliLearn(il));
    end
end