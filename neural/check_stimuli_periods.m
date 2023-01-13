function bad_period = check_stimuli_periods(periods)
% Check if the length of image durations are correct
    bad_period = 0;
    if range(periods)>=20000
        fprintf('%s\n',['Periods are  variable !!! ' num2str(range(periods))]);
        bad_period=1;
    end
    if mean(periods)>1100*1000
        fprintf('This session used a longer learning stimulus duration %d\n',...
            round(mean(periods)/1000));
    end
end

