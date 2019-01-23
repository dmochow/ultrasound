clear all; close all; clc

%%
rStr='r010';
freqStr='10Hz';
condStr='Active';
fs=30000;
dsr=5; %downsampling ratio
flHz=300; % lowest frequency that we're interested in
fhHz=3000; % highest frequency we interested in
butterOrder=3;
probeMap=[17,16,18,15,19,14,20,13,...
    21,12,22,11,23,10,24,9,...
    25,8,26,7,27,6,28,5,...
    29,4,30,3,31,2,32,1]';

dataFolder=['../data/' rStr '/*.mat' ];
filenames=dir(dataFolder);
nFiles=numel(filenames);
for f=1:nFiles
    tFilename=filenames(f).name;
    if contains(tFilename,'Pre') && contains(tFilename,freqStr) && contains(tFilename,condStr)
        preFilename=tFilename;
    elseif contains(tFilename,'Dur') && contains(tFilename,freqStr) && contains(tFilename,condStr)
        periFilename=tFilename;
    elseif contains(tFilename,'Pos') && contains(tFilename,freqStr) && contains(tFilename,condStr)
        postFilename=tFilename;
    else
    end
end

%%
matFilenames{1}=fullfile(dataFolder,preFilename);
matFilenames{2}=fullfile(dataFolder,periFilename);
matFilenames{3}=fullfile(dataFolder,postFilename);
samples={60*0*fs+1:60*3*fs;60*0*fs+1:60*3*fs;60*0*fs+1:60*3*fs;60*0*fs+1:60*3*fs;60*0*fs+1:60*3*fs;60*0*fs+1:60*3*fs}; % 3 min 

%%
nFilenames=numel(matFilenames);
for f=1:nFilenames
    f
    load(matFilenames{f}); 
    data(:,:,f)=amplifier_data(:,samples{f}); clear amplifier_data;
end

%%
fl=flHz/(fs/2);
fh=fhHz/(fs/2);
[b,a]=butter(butterOrder,[fl fh]);
data=filter(b,a,data,[],2);
nHarm=floor(fhHz/60); 
for h=1:nHarm
    h
    [bn,an]=butter(2,[h*60-2 h*60+2]/fs*2,'stop');
    data=filter(bn,an,data,[],2);
end


% Hn=dfilt.df2t(bn,an);
% Hf=dfilt.cascade(Hf,Hn);
% data=filter(Hf,data,2);

%%
% downsample data (assumes lpf already done)
dataDown=data(:,1:dsr:end,:);
fsr=fs/dsr;

%% reject gross artifacts
Q=40; % how many standard deviations above median(abs)/0.6745 to call it artifact
xtent=0.1*fsr; % how many adjacent samples to remove alongside artifacts
for e=1:size(dataDown,3)
    e
    dataDownClean(:,:,e) = zeroBadSamples(dataDown(:,:,e),Q,xtent);
end

%% spike detection
threshStds=5; fs=fsr; type='neg';
% [spikesPreSh,waveformPreSh] = detectSpikes(dataDownClean(:,:,1),threshStds,fs,type);
% [spikesPeriSh,waveformPeriSh] = detectSpikes(dataDownClean(:,:,2),threshStds,fs,type);
% [spikesPostSh,waveformPostSh] = detectSpikes(dataDownClean(:,:,3),threshStds,fs,type);
[spikesPre,waveformPre,muaPre] = detectSpikes(dataDownClean(:,:,1),threshStds,fs,type);
[spikesPeri,waveformPeri,muaPeri] = detectSpikes(dataDownClean(:,:,2),threshStds,fs,type);
[spikesPost,waveformPost,muaPost] = detectSpikes(dataDownClean(:,:,3),threshStds,fs,type);

%%
% count
% [nSpikes(:,1),b]=cellfun(@size,spikesPreSh);
% [nSpikes(:,2),b]=cellfun(@size,spikesPeriSh);
% [nSpikes(:,3),b]=cellfun(@size,spikesPostSh);
[nSpikes(:,1),~]=cellfun(@size,spikesPre);
[nSpikes(:,2),~]=cellfun(@size,spikesPeri);
[nSpikes(:,3),~]=cellfun(@size,spikesPost);
spikeRate=nSpikes/180; % TODO: generalize duration

% %%
chanInds=[26 25 8 9 7 17 16 2 1 32];
nRows=2; nCols=5;
figure;
for ch=1:numel(chanInds)
    hs(ch)=subplot(nRows,nCols,ch);
    plot(spikeRate(chanInds(ch),:),'--o');
end

% %%
% delSpikeRateSh=spikeRate(:,3)-spikeRate(:,1);
% delSpikeRate=spikeRate(:,6)-spikeRate(:,4);
figure; subplot(121);
imagesc(delSpikeRateSh(probeMap));
colormap('jet');
% caxis([-5 5]);
% subplot(122);
% imagesc(delSpikeRate(probeMap));
% colormap('jet');
% caxis([-5 5]);
% 
% %%
% targetInds=chanInds(1:5);
% nTargetInds=setdiff(1:32,targetInds);
% figure
% subplot(221);
% bar([1 2],[ mean(spikeRate(targetInds,4)) mean(spikeRate(targetInds,6))]);
% set(gca,'XTickLabel',{'Pre','Post'});
% title('Target Active');
% subplot(222);
% bar([1 2],[ mean(spikeRate(nTargetInds,4)) mean(spikeRate(nTargetInds,6))]);
% set(gca,'XTickLabel',{'Pre','Post'});
% title('Non-Target Active');
% subplot(223);
% bar([1 2],[ mean(spikeRate(targetInds,1)) mean(spikeRate(targetInds,3))]);
% set(gca,'XTickLabel',{'Pre','Post'});
% title('Target Sham');
% subplot(224);
% bar([1 2],[ mean(spikeRate(nTargetInds,1)) mean(spikeRate(nTargetInds,3))]);
% set(gca,'XTickLabel',{'Pre','Post'});
% title('Non-Target Sham');
% print -dpng ../figures/spikeRate_prePeriPost_40Hz

%%
save(['../data/precomputed/' rStr '_prePeriPost_spikes_' condStr '_' freqStr],...
    'spikesPre','waveformPre', ...
    'spikesPeri','waveformPeri', ...
    'spikesPost','waveformPost', ...
    'muaPre','muaPeri','muaPost'); 

%% smooth firing rate
%%
% %% now apply sliding window to get firing rate
% % kernel for firing rate computation
% sigma=0.01; % for computing firing rate
% nKernel=30*70;
% tKernel=(-nKernel/2+1:1:nKernel/2)/fsr;
% kernel=1/sqrt(2*pi*sigma.^2)*exp(-(tKernel.^2/(2*sigma^2)));
% nPad=sum(tKernel<0);
% figure; plot(tKernel,kernel,'k');
% 
% %%
% fr=zeros(size(mua));
% for ch=1:size(mua,1)
%     for bl=1:size(mua,3)
%         paddedData=[zeros(1,nPad) mua(ch,:,bl)];
%         smoothedData=conv(paddedData,kernel,'same');
%         fr(ch,:,bl)=smoothedData(nPad+1:end);
%     end
% end

