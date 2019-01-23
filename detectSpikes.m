function [spikes,waveform,mua,thresh] = detectSpikes(data,threshStds,fs,type,thresh)
%detect spikes based on amplitude thresholding
% 

if nargin<5, thresh=''; end
if nargin<4, type='neg'; end
if nargin<3, fs=30000; end
nChannels=size(data,1);
nSamples=size(data,2);
chMeds=median(abs(data),2);

if isempty(thresh)
    thresh=chMeds/0.6745*threshStds;
end
%threshArt=chMeds/0.6745*25; %25 standard deviations

switch type
    case 'pos'
        mdiff=cat(2,zeros(size(data,1),1),diff(data>repmat(thresh,1,nSamples),1,2));
        %madiff=cat(2,zeros(size(data,1),1),diff(data>repmat(threshArt,1,nSamples),1,2));
    case 'neg'
        mdiff=cat(2,zeros(size(data,1),1),diff(data<repmat(-thresh,1,nSamples),1,2));
        %madiff=cat(2,zeros(size(data,1),1),diff(data<repmat(-threshArt,1,nSamples),1,2));
end

%
[ch,tm]=find(mdiff>0); % get channel and sample
spikeTable=[ch, tm];

%
spikeInds=mdiff>0;
mua=zeros(size(data));
mua(spikeInds)=1;  


%artInds=find(madiff>0);
%mua(artInds)=0;

% now compute and store waveforms
Tspike=round(fs/500); % 2 ms
waveform=cell(nChannels,1);
spikes=cell(nChannels,1);

% for ch=1:nChannels
%     spikeTimes=find(mua(ch,:));
%     for s=1:numel(spikeTimes)
%         if spikeTimes(s)+Tspike/2 <= nSamples
%             %waveform{ch}(:,s)=data(ch,spikeTimes(s)+1:spikeTimes(s)+round(fs/500)); % store 2 ms
%             waveform{ch}(:,s)=data(ch,spikeTimes(s)-Tspike/2+1:spikeTimes(s)+Tspike/2); % store 2 ms
%         end
%     end
% end

% spikes=[];
% waveform=[];

for s=1:size(spikeTable,1)
    %s/size(spikeTable,1)
    chIdx=spikeTable(s,1);
    tmSample=spikeTable(s,2);
    if tmSample+Tspike/2 <= nSamples & tmSample-Tspike/2+1 > 0
        waveform{chIdx}=cat(1,waveform{chIdx},data(chIdx,tmSample-Tspike/2+1:tmSample+Tspike/2)); % store 2 ms
    end
    
    spikes{chIdx}=cat(1,spikes{chIdx},tmSample);
end




end

