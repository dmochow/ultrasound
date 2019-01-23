clear all; close all; clc


filenames={'../data/precomputed/r008_prePeriPost_spikes_sham10Hz',...
    '../data/precomputed/r008_prePeriPost_spikes_active10Hz',...
    '../data/precomputed/r008_prePeriPost_spikes_sham40Hz',...
    '../data/precomputed/r008_prePeriPost_spikes_active40Hz'};
DURATION_SECS=180;
targetChanInds=[26 25 8 7 9];
nonTargetChanInds=[17 16 1 2 32];
%nonTargetChanInds=setdiff(1:32,targetChanInds);
nFilenames=numel(filenames);
titleStrings={'Sham 10 Hz','Active 10 Hz','Sham 40 Hz','Active 40 Hz'};
probeMap=[17,16,18,15,19,14,20,13,...
    21,12,22,11,23,10,24,9,...
    25,8,26,7,27,6,28,5,...
    29,4,30,3,31,2,32,1]';
probeMap=probeMap(end:-1:1);


%% target contacts
figure(1)
for f=1:nFilenames
    load(filenames{f});
    [nSpikes(:,1),~]=cellfun(@size,spikesPre);
    [nSpikes(:,2),~]=cellfun(@size,spikesPeri);
    [nSpikes(:,3),~]=cellfun(@size,spikesPost);
    spikeRate=nSpikes/DURATION_SECS;
    allSpikeRate(:,:,f)=spikeRate;
    % firing rate
%     fr{f,1} = spikeTrain2fr(muaPre);
%     fr{f,2} = spikeTrain2fr(muaPeri);
%     fr{f,3} = spikeTrain2fr(muaPost);
    
    barHeights=mean(spikeRate(targetChanInds,[1 3]));
    hs(f)=subplot(2,2,f);
    bar([1 2],barHeights);
    title(titleStrings{f});
    set(gca,'Xticklabel',{'Pre','Post'});
end

%% non-target contacts
figure(2)
for f=1:nFilenames
    load(filenames{f});
    [nSpikes(:,1),~]=cellfun(@size,spikesPre);
    [nSpikes(:,2),~]=cellfun(@size,spikesPeri);
    [nSpikes(:,3),~]=cellfun(@size,spikesPost);
    spikeRate=nSpikes/DURATION_SECS;
    barHeights=mean(spikeRate(nonTargetChanInds,[1 3]));
    hs(f)=subplot(2,2,f);
    bar([1 2],barHeights);
    title(titleStrings{f});
    set(gca,'Xticklabel',{'Pre','Post'});
end

%%
figure; hold on ;
for f=1:nFilenames
subplot(2,2,f);
barh(1:32,allSpikeRate(probeMap,3,f)-allSpikeRate(probeMap,1,f));
end
