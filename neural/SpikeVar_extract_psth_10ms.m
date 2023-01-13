
clear;load ProcessedOutput ; %this is the output of Noneural_main_release.m
od=pwd; SwitchFolders('OurCode');

%from table 1 of data descriptor - ERROR in original, since session id
%order in table 1 is not in sequence
addpath(genpath('/Users/waschke/OneDrive/Old data and code/Old data and code/RecogMemory_MTL_release_v2/'))
[NOsessions, NO_listOf_allUsable] = defineNOsessions_release(); %use NOind to loopinto this
agevec=defineAges;
cval=NO_listOf_allUsable; % redundant
stimdir = '/Users/waschke/OneDrive/Old data and code/SV/RecogMemory_MTL_release_v2/Code/dataRelease/stimFiles/';
%%

%output of findComplexities.m
load([stimdir 'SpikeVar_stiminfo_V1.mat']);
load([stimdir 'SpikeVar_stiminfo_V2']);
load([stimdir 'SpikeVar_stiminfo_V3']);
Lstim=[v1_stims.learnStim;v2_stims.learnStim;v3_stims.learnStim];
Lstim2=[v1_stims.learnStimRepeated;v2_stims.learnStimRepeated];
Lcat=[v1_stims.learnStimCat;v2_stims.learnStimCat;v3_stims.learnStimCat];
Lcat2=[v1_stims.learnStimRepeatedCat;v2_stims.learnStimRepeatedCat];
Lcatnum=[v1_stims.learnStimCatNum;v2_stims.learnStimCatNum;v3_stims.learnStimCatNum];
Lcat2num=[v1_stims.learnStimRepeatedCatNum;v2_stims.learnStimRepeatedCatNum];

NOsessions=NOsessions(cval);
varval=[NOsessions.variant];

MatStruct = struct('participant',NaN,'session',NaN,'taskvariant',NaN,...
    'NOind',NaN,'numtrials',NaN,'diagnosisCode',NaN,...
    'trial',NaN,'stimuliLearn',NaN,'LearnCorrect',NaN,'LearnRT',NaN,...
    'RecogResponse',NaN,'RecogRT',NaN,...
    'neuronData',[],'age',NaN,'LearnImageDur',NaN,...
    'stimuliLearnCat',NaN,'stimuliLearnCatNum',NaN);


MatStruct =repmat(MatStruct,1,length(allData));
ltl=NaN(1,length(allData));%just to use later in dropme
badperiods=ltl;
NeuronL=1;
%%
for ind=1:length(allData)
    fprintf('%d\n',ind);
    %**********SESSION LEVEL************
    %in next step, repmat these for matrix, the concatenate
    MatStruct(ind).participant=NOsessions(ind).patient_nb;
    MatStruct(ind).session=NOsessions(ind).patientsession;
    MatStruct(ind).taskvariant=NOsessions(ind).variant;
    MatStruct(ind).diagnosisCode=NOsessions(ind).diagnosisCode;
    MatStruct(ind).numtrials=length(allData(ind).stimuliLearn); %numtrials.. just for safety
    MatStruct(ind).NOind=cval(ind); %index into defineNoSessions.. just for safety
    MatStruct(ind).age=agevec(ind);
    %**********Across NEURONS, TRIAL LEVEL*****************
    %vector form
    MatStruct(ind).trial=[1:length(allData(ind).stimuliLearn)]'; %sometimes this is 50
    MatStruct(ind).stimuliLearn=allData(ind).stimuliLearn'; %this is actually not unique..changing below
    ltl(ind)=length(allData(ind).LresponsesEventfile); %just to use later in dropme
    MatStruct(ind).LearnCorrect=allData(ind).learnAccurate;
    MatStruct(ind).LearnRT=allData(ind).LRTRecog_raw;
    
    %do a check of stimulus order here
    if length(allData(ind).stimuliLearn)==100
        if any(allData(ind).stimuliLearn~=Lstim(varval(ind),:)),keyboard;end
        stimuliLearnCat=Lcat(varval(ind),:)';
        stimuliLearnCatNum=Lcatnum(varval(ind),:)';
    elseif length(allData(ind).stimuliLearn)==50
        if any(allData(ind).stimuliLearn~=Lstim2(varval(ind),:)),keyboard;end
        stimuliLearnCat=Lcat2(varval(ind),:)';
        stimuliLearnCatNum=Lcat2num(varval(ind),:)';
    end
    
    MatStruct(ind).stimuliLearnCat=stimuliLearnCat;
    MatStruct(ind).stimuliLearnCatNum=stimuliLearnCatNum;
    
    if stimuliLearnCatNum>99, keyboard;end
    MatStruct(ind).stimuliLearn=(stimuliLearnCat+10)*100+stimuliLearnCatNum; %4 digit stimulus ID, first two digits are category + 10
    
    %Now finding the recognition response and recognition RT
    % for the repeated learning stimuli
    sl=allData(ind).stimuliLearn;
    sr=allData(ind).stimuliRecog;
    rr=allData(ind).recogResponses;
    rt=allData(ind).RTRecog_raw;
    
    RecogResponse=NaN(1,length(sl));
    RecogRT=RecogResponse;
    [~,repeatind]=ismember(sl,sr);
    RecogResponse(repeatind~=0)=rr(repeatind(repeatind~=0));
    RecogRT(repeatind~=0)=rt(repeatind(repeatind~=0));
    MatStruct(ind).RecogResponse=RecogResponse';
    MatStruct(ind).RecogRT=RecogRT';
    
    curLearnImageOn=allData(ind).periodsLearn0(:,2:3);%starts with image onset and ends at image offset
    if range(diff(curLearnImageOn,1,2))>=20000,fprintf('%s\n',['Periods are  variable !!! '...
            num2str(range(diff(curLearnImageOn,1,2)))]);badperiods(ind)=1;
    end
    if mean(diff(curLearnImageOn,1,2))>1100*1000,fprintf(...
            'This session used a longer learning stimulus duration %d\n',...
            round(mean(diff(curLearnImageOn,1,2))/1000));
    end
    
    if mean(diff(curLearnImageOn,1,2))<1100*1000,MatStruct(ind).LearnImageDur=1000;
    else MatStruct(ind).LearnImageDur=2000;
    end
    
    
    %********************NOW NEURON * TRIAL LEVEL*************************
    curspikes=allData(ind).cellStats;
    curtotstats=allData(ind).allStats;
    spikeMat=struct('Neuron',NaN,'brainArea',NaN);%one for each neuron
    %'MU',NaN
    spikeMat=repmat(spikeMat,1,length(curspikes));
    
    for spikeind=1:length(curspikes)
        numTrials=size(allData(ind).periodsLearn,1);
        spikeMat(spikeind).Neuron=repmat(spikeind,numTrials,1);%neuron number
        spikeMat(spikeind).NeuronL=repmat(NeuronL,numTrials,1);NeuronL=NeuronL+1; %running neuron count
        spikeMat(spikeind).brainArea=repmat(curspikes(spikeind).brainAreaOfCell,numTrials,1);
        spikeMat(spikeind).trialNumber=(1:numTrials)'; %should always be 100 or 50, usually 100.
        curSpikeTrain=curspikes(spikeind).timestamps;
        
        %Use findCounts.m to count spikes in the spiketrain
        spikeMat(spikeind).clusterID=repmat(curspikes(spikeind).origClusterID,numTrials,1);
        %         a test for an on-response to the learning image
        curbinwidth=1000; %milliseconds
        precurbinwidth=1000;postcurbinwidth=1500;
        basecounts=findCounts(curSpikeTrain,[curLearnImageOn(:,1)-(precurbinwidth*1000)...
            curLearnImageOn(:,1)],precurbinwidth*1000);% multiplier to get away from microseconds
        respcounts=findCounts(curSpikeTrain,[curLearnImageOn(:,1)...
            curLearnImageOn(:,1)+(postcurbinwidth*1000)],postcurbinwidth*1000);
        basecounts=basecounts/(precurbinwidth/1000);
        respcounts=respcounts/(postcurbinwidth/1000);
        spikeMat(spikeind).basecounts=basecounts;
        spikeMat(spikeind).respcounts=respcounts;
        spikeMat(spikeind).IncreaseStatus=repmat(nanmean(respcounts)>nanmean(basecounts),length(basecounts),1);
        
        pval=pairedPermute(basecounts,respcounts,1000,0);
        spikeMat(spikeind).ONpval_perm=repmat(pval,numTrials,1);
        [dum,pval]=ttest(basecounts,respcounts);
        spikeMat(spikeind).ONpval=repmat(pval,numTrials,1);
        pval=signrank(basecounts,respcounts);
        spikeMat(spikeind).sONpval=repmat(pval,numTrials,1);
        
        pval=pairedPermute(basecounts,respcounts,1000,1);
        spikeMat(spikeind).ONpval_perm2=repmat(pval,numTrials,1);
        [dum,pval]=ttest(basecounts,respcounts,'tail','left');
        spikeMat(spikeind).ONpval2=repmat(pval,numTrials,1);
        
        % neuronal stats from the recognition phase, based on original
        % release
        spikeMat(spikeind).rONpval=repmat(curtotstats(spikeind,4),numTrials,1); 
        spikeMat(spikeind).rMS=repmat(curtotstats(spikeind,5),numTrials,1);
        spikeMat(spikeind).rVS=repmat(curtotstats(spikeind,6),numTrials,1);
        spikeMat(spikeind).rONpval0=repmat(curspikes(spikeind).pVals0,numTrials,1);
        
        % Use 10 ms bins to count spikes
        curbinwidth=10;
        ocurLearnImageOn=curLearnImageOn;
        ocurLearnImageOn(:,2)=curLearnImageOn(:,2)+500000; 
        
        psth=findCounts(curSpikeTrain,ocurLearnImageOn,curbinwidth*1000);
        opsth=psth;%this is the full-sized psth
        cut=round(500000/(curbinwidth*1000));
        psth=psth(:,1:(end-cut)); %removing last x bins
        
        spikeMat(spikeind).mean10=nanmean(psth,2);
        
        numcol=round(2500/curbinwidth);numrows=size(opsth,1);numpsthcol=size(opsth,2);
        spikeMat(spikeind).psth10=[opsth NaN(numrows,numcol-numpsthcol)];
        spikeMat(spikeind).psth10(:,round(2500/curbinwidth):end) = [];
        
        ocurLearnImageOn(:,1)=curLearnImageOn(:,1)-(1000*1000);
        opsth=findCounts(curSpikeTrain,ocurLearnImageOn,curbinwidth*1000);%the multiplier is because in microseconds
        numcol=round(3500/curbinwidth);numrows=size(opsth,1);numpsthcol=size(opsth,2);
        spikeMat(spikeind).opsth10=[opsth NaN(numrows,numcol-numpsthcol)];
        spikeMat(spikeind).opsth10(:,round(3500/curbinwidth):end) = [];
      
    end
    MatStruct(ind).neuronData=spikeMat;
end
%%
fprintf('%s\n','');

%e drop participants with only 89 learning trials
dropme=([NOsessions.blockIDRecog]==3|ltl==89); %|badperiods==1%|badperiods==1
MatStruct=MatStruct(~dropme);
NOsessions=NOsessions(~dropme);

totnum=0;
for ind=1:length(MatStruct)
    for ind1=1:length(MatStruct(ind).neuronData)
        totnum=totnum+1;
    end
end
onpval=NaN(1,totnum);
onpval2=onpval;
ronpval=onpval;
onpval_perm=onpval;
onpval_perm2=onpval;
sonpval=onpval;
ronpval0=onpval;
totnum=1;
for ind=1:length(MatStruct)
    for ind1=1:length(MatStruct(ind).neuronData)
        onpval(totnum)=MatStruct(ind).neuronData(ind1).ONpval(1);
        onpval2(totnum)=MatStruct(ind).neuronData(ind1).ONpval2(1);
        ronpval(totnum)=MatStruct(ind).neuronData(ind1).rONpval(1);
        ronpval0(totnum)=MatStruct(ind).neuronData(ind1).rONpval0(1);
        onpval_perm(totnum)=MatStruct(ind).neuronData(ind1).ONpval_perm(1);
        onpval_perm2(totnum)=MatStruct(ind).neuronData(ind1).ONpval_perm2(1);
        sonpval(totnum)=MatStruct(ind).neuronData(ind1).sONpval(1);
        totnum=totnum+1;
    end
end

%Here construct the matrix
fulltable=table();
for ind=1:length(MatStruct)
    curmat=MatStruct(ind);
    curspikes=MatStruct(ind).neuronData;
    
    spiketable=table();
    for neuronind=1:length(curspikes)
        spiketable=[spiketable;struct2table(curspikes(neuronind))];
    end
    slen=size(spiketable,1);
    otable=struct2table(struct('Participant',repmat(MatStruct(ind).participant,1,slen)',...
        'Session',repmat(MatStruct(ind).session,1,slen)','TaskVariant',...
        repmat(MatStruct(ind).taskvariant,1,slen)','NumTrials',...
        repmat(MatStruct(ind).numtrials,1,slen)','DiagnosisCode',...
        repmat(MatStruct(ind).diagnosisCode,1,slen)','Age',...
        repmat(MatStruct(ind).age,1,slen)','learnImageDur',...
        repmat(MatStruct(ind).LearnImageDur,1,slen)'));
    
    otable2=struct2table(struct('stimuliLearn',...
        repmat(MatStruct(ind).stimuliLearn,length(curspikes),1),...
        'stimuliLearnCat',repmat(MatStruct(ind).stimuliLearnCat,...
        length(curspikes),1),'stimuliLearnCatNum',...
        repmat(MatStruct(ind).stimuliLearnCatNum,length(curspikes),1),...
        'learnCorrect',repmat(MatStruct(ind).LearnCorrect,length(curspikes),1),...
        'learnRT',repmat(MatStruct(ind).LearnRT,length(curspikes),1),...
        'RecogResponse',repmat(MatStruct(ind).RecogResponse,length(curspikes),1),...
        'RecogRT',repmat(MatStruct(ind).RecogRT,length(curspikes),1)));
    
    fulltable=[fulltable; [otable otable2 spiketable]];
end
%%
save SPSSinput_10_msec fulltable MatStruct
cd(od);
writetable(fulltable, 'SPSSinput_10_msec.csv');






