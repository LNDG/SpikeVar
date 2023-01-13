function get_meta_data()
    [NOsessions, NO_listOf_allUsable] = defineNOsessions_release();
    for c = 1:120
        NOsessions(c).index = c;
    end
    % delete not useable sessions
    NOsessions = NOsessions(NO_listOf_allUsable);
    base_dir = '/Users/kamp/PhD/spikevar';
    file_name = fullfile(base_dir,'data', 'neuro', 'NOsessions.mat');
    save(file_name, 'NOsessions')
end