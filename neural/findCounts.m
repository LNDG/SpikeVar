function [outmat,ax]=findCounts(timeStamps,periods,binwidth,overlap)
%function outmat=findCounts(timeStamps,periods,binwidth)
%size(periods,1) is the number of trials
%spikes are counted left-continuous.. so >=left bin and <right bin, based
%on histc behavior
%if non-overlapping (last value is 0 or NaN or does not exist), histc will be used for
%speed.
%can also use for counts by setting binwidth equal to total count-duration
%HUGE assumption that periods is increasing..have to handle backwards time
%before function call

if ~exist('overlap','var')||overlap<=0||isnan(overlap)
    overlap=0;
else overlap=0; fprintf('%s\n','Sorry.. non-overlapping bins not handled yet');
end


% if range(diff(periods,1,2))>20000,badrange=1;else badrange=0;end %periods vary by more than 20 ms.
%***watch out ! diff has weird input arguments
periodDur=nanmean(diff(periods,1,2)); %taking mean period duration... periods of different lengths ignored
periodDur=round(periodDur/binwidth)*binwidth; %rounding to nearest multiple of binwidth, needed due to microseconds
ax=0:binwidth:periodDur;%will factor in correctly

outmat=NaN(size(periods,1),length(ax)-1);
for ind=1:size(periods,1)
    temp=histc(timeStamps,periods(ind,1)+ax);
    temp=temp(1:(end-1)); %ignoring last bin
    outmat(ind,:)=temp;
end