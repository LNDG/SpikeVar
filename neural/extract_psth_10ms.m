clear
clc
% check if we are in the correct directory, change it if needed
base_dir = '/Users/kamp/PhD/spikevar';
neuro_dir = fullfile(base_dir, 'repo', 'SpikeVar', 'neural/');
data_dir = fullfile(base_dir,'data', 'neuro/');
stim_dir = fullfile(base_dir,'stimuli_rn/');

load([stim_dir 'SpikeVar_stiminfo_V1.mat'], 'v1_stims');
load([stim_dir 'SpikeVar_stiminfo_V2.mat'], 'v2_stims');
load([stim_dir 'SpikeVar_stiminfo_V3.mat'], 'v3_stims');
stim_info = struct('v1', v1_stims, 'v2', v2_stims, 'v3', v3_stims);

load([data_dir 'ProcessedOutput.mat'], 'allData');
load([data_dir 'NOsessions.mat'], 'NOsessions');

% check if sessionIDs are the same in allData and NOsessions
if ~all(cellfun(@strcmp, {allData.sessionID}, {NOsessions.sessionID}))
    keyboard
end

age = get_ages();

%%
% initiat data structure
num_data = length(allData); 
neuron_count = 1; %running neuron count

spike_data = struct(...
    'participant',NaN,'session',NaN,'taskvariant',NaN,...
    'NOind',NaN,'numtrials',NaN,'diagnosisCode',NaN,...
    'trial',NaN,'stimuliLearn',NaN,'LearnCorrect',NaN,'LearnRT',NaN,...
    'RecogResponse',NaN,'RecogRT',NaN,...
    'neuronData',[],'age',NaN,'LearnImageDur',NaN,...
    'stimuliLearnCat',NaN,'stimuliLearnCatNum',NaN);

for i=1:num_data
    %**********SESSION LEVEL************
    spike_data(i).participant = NOsessions(i).patient_nb;
    spike_data(i).session = NOsessions(i).patientsession;
    task_variant = NOsessions(i).variant;
    spike_data(i).task_variant = task_variant;
    spike_data(i).diagnosisCode = NOsessions(i).diagnosisCode;
    spike_data(i).NOind = NOsessions(i).index; %index in NoSessions
    num_trials = length(allData(i).stimuliLearn);
    spike_data(i).num_trials = num_trials;
    spike_data(i).age = age(i, 2);
    spike_data(i).bad_periods = NaN;
    spike_data(i).LresponsesEventfile = length(allData(i).LresponsesEventfile); % used to drop later

    
    %**********Across NEURONS, TRIAL LEVEL*****************
    spike_data(i).trial = [1:num_trials]'; % num_trials is 100 or 50
    spike_data(i).stimuliLearn = allData(i).stimuliLearn'; %this is actually not unique..changing below
    spike_data(i).LearnCorrect = allData(i).learnAccurate;
    spike_data(i).LearnRT = allData(i).LRTRecog_raw;
    
    % create unique stimulus codes
    if num_trials == 100
        % check of stimulus order
        if any(spike_data(i).stimuliLearn ~= stim_info.(sprintf('v%d', task_variant)).learnStim')
            keyboard
        end
        stimuliLearnCat = stim_info.(sprintf('v%d', task_variant)).learnStimCat';
        stimuliLearnCatNum = stim_info.(sprintf('v%d', task_variant)).learnStimCatNum';
    elseif num_trials == 50
        if any(spike_data(i).stimuliLearn ~= stim_info.(sprintf('v%d', task_variant)).learnStimRepeated')
            keyboard
        end  
        stimuliLearnCat = stim_info.(sprintf('v%d', task_variant)).learnStimRepeatedCat';
        stimuliLearnCatNum = stim_info.(sprintf('v%d', task_variant)).learnStimRepeatedCatNum';
    end
    
    if any(stimuliLearnCat > 99)
        keyboard
    end
    
    spike_data(i).stimuliLearnCat = stimuliLearnCat;
    spike_data(i).stimuliLearnCatNum = stimuliLearnCatNum;
    
    %4 digit stimulus code, first two digits are category + 10
    spike_data(i).stimuliCode =(stimuliLearnCat + 10) * 100 + stimuliLearnCatNum;
    
    % Find the recognition response and recognition RT
    % for the repeated learning stimuli
    RecogResponse = NaN(1, num_trials);
    RecogRT = NaN(1, num_trials);
    [~,repeated_ind] = ismember(allData(i).stimuliLearn, allData(i).stimuliRecog);
    RecogResponse(repeated_ind ~= 0)= allData(i).recogResponses(repeated_ind(repeated_ind~=0));
    RecogRT(repeated_ind~=0)=allData(i).RTRecog_raw(repeated_ind(repeated_ind~=0));
    
    spike_data(i).RecogResponse=RecogResponse';
    spike_data(i).RecogRT=RecogRT';
    
    % image on and off 
    LearnImageOnOff = allData(i).periodsLearn0(:,2:3);
    % time periods that each stimulus is shown
    periods = diff(LearnImageOnOff,1,2);  
    spike_data(i).bad_periods = check_stimuli_periods(periods);
    if mean(periods)<1100*1000
        spike_data(i).LearnImageDur=1000;
    else
        spike_data(i).LearnImageDur=2000;
    end
    
    %******************** NEURON * TRIAL LEVEL*************************
    cell_stats = allData(i).cellStats;
    n_cells = length(cell_stats);
    all_stats = allData(i).allStats;
    cell_struct = struct('Neuron',NaN,'brainArea',NaN);
    
    for c_idx = 1:length(cell_stats)
        cell_struct(c_idx).Neuron = c_idx;
        cell_struct(c_idx).NeuronCount = neuron_count;
        neuron_count = neuron_count+1; 
        cell_struct(c_idx).brainArea= cell_stats(c_idx).brainAreaOfCell;
        cell_struct(c_idx).trialNumber = (1:num_trials)'; 
        cell_struct(c_idx).clusterID = cell_stats(c_idx).origClusterID;
        
        % get spike train
        spike_train = cell_stats(c_idx).timestamps;
               
        % add 500ms to the end of each period
        LearnImageOnOff(:,2)=LearnImageOnOff(:,2)+ 500*1000; 
        % get binned spike counts for every stimulus with bin width of 10ms
        bin_width = 10; 
        stimuliSpikeCounts = findCounts(spike_train, LearnImageOnOff, bin_width*1000);
        
        % pad table with NaN
        num_col =round(2500/bin_width);
        num_rows=size(stimuliSpikeCounts,1);
        cell_struct(c_idx).stimuliSpikeCounts=[stimuliSpikeCounts NaN(num_rows, num_col-size(stimuliSpikeCounts,2))];
        cell_struct(c_idx).stimuliSpikeCounts(:,round(2500/bin_width):end) = [];
        
        %removing last x bins
        cut = round(500000/(bin_width*1000));
        stimuliSpikeCounts = stimuliSpikeCounts(:,1:(end-cut)); 
        cell_struct(c_idx).meanStimuliSpikeCounts = nanmean(stimuliSpikeCounts,2);        
    end
    spike_data(i).neuronData=cell_struct;
end

% drop participants with only 89 learning trials
drop_me =([NOsessions.blockIDRecog]==3|[spike_data.LresponsesEventfile]==89); %|badperiods==1%|badperiods==1
spike_data=spike_data(~drop_me);

% save data struct
save([data_dir 'SpikeVar_spike_data.mat'], 'spike_data')



