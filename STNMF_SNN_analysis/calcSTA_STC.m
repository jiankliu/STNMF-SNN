function [STA, STC] = calcSTA_STC(stim, train, maxLag, sampFreq)
%	CALCULATE STA AND STC
%
% [STA, STC] = calcSTA(stim, train, maxLag)
%
%ARGS:
% stim - stimulus
% train - binary spike train 
%         (i.e. 0 if no spike and 1 if spike at [t,t+dt]) 
%         and having the same rate as the stimulus
% maxLag - maximum lag in pts
%
%RETURNS:
% BOTH VALUES ARE NOT NORMALIZED
% STA - spike triggered average of size [1*maxLag+1,1]
% STC - spike triggered covariance matrix of size [1*maxLag+1,1*maxLag+1],
%
%NOTE:
% stim and train must have the same sampling frequency
% STC is only an approximation as i'm calc'ing it in chunks of 1000 spikes
% (!!!!EXPERIMENTAL!!!! plus needs much more data than STA due to higher dimensionality)
%TODO:
% allow variable sampling rates
%
% created 05/08/09, JC
% lastmod 06/06/28, JC

spikes = find(train==1)
spikes = spikes(find(spikes>maxLag*sampFreq/1000 & spikes<length(train)-maxLag*sampFreq/1000));

%_______ S T A and S T C _______________________________________
fprintf('\ndoing STA and STC\t');
slice = 1:maxLag*sampFreq/1000;
dt = 1;
STA = zeros(1,length(slice));
STC = zeros(floor(length(slice)/dt))

%___ online mean calc _____________________________________
chunk = 50;
idx = [0,chunk*(1:ceil(length(spikes)/chunk)-1), length(spikes)];
for ch = 2:length(idx)
    range = idx(ch-1)+1:idx(ch);
    trigIdx = repmat(slice,length(range),1)+repmat(spikes(range),1,length(slice));
    trigStim = stim(trigIdx);
    m = mean(trigStim);
    STA = STA + m;
end

STA = STA./length(idx);%normalize STA by spike count

%___ online covariance calc _____________________________________
% subSlice = dt:dt:length(STA);
% for ch = 2:length(idx)
%     fprintf('.');
%     range = idx(ch-1)+1:idx(ch);
%     trigIdx = repmat(subSlice,length(range),1)+repmat(spikes(range),1,length(subSlice));
%     trigStim = stim(trigIdx);
%     trigStim = trigStim - repmat(STA(subSlice),size(trigStim,1),1);
%     STC = STC + trigStim'*trigStim;
% end
% 
% STC = STC./length(idx);%normalize STC - correctly???
fprintf('\n');