clear all; close all; clc
% this script assumes that rhd file has been read in by 
%read_Intan_RHD2000_file.m and then save to a mat file which 
%serves as the starting point of this script
%
% (c) Jacek "Loquacious D" Dmochowski, 2018-

matFilename='../data/r003/smartboxFirstData.mat';

%% 
% grab the data and basic parammies
load(matFilename);
data=amplifier_data;
nChannels=size(amplifier_data,1);
nSamples=size(amplifier_data,2);
fs=frequency_parameters.amplifier_sample_rate;

% specify some other 'rammies
channelsKeep=1:nChannels; % subset of channels that we want to keep
flHz=300; % lowest frequency that we're interested in
fhHz=4000; % highest frequency we interested in
fsr=floor(fs/fhHz);
fsd=fs/fsr;
butterOrder=4;
dimFilter=2; % time is usually the second (column) dimension
fl=flHz/(fs/2); % trust me
fh=fhHz/(fs/2);
[b,a]=butter(butterOrder,[fl fh]);
data=filter(b,a,data,[],2);

% notch filter to remove 60 Hz + harmonics
nHarm=50;
for h=1:nHarm
    [bn,an]=butter(2,[h*60-2 h*60+2]/fs*2,'stop');
    data=filter(bn,an,data,[],2);
end

%data=downsample(data.',fsr).';

%% 
% spike detection parameters
threshStds=5; % spike detection threshold (multiple of median(abs)/0.6745)
[spikes,waveform,mua,thresh] = detectSpikes(data(:,:),threshStds,fs,'neg');

%%
% pick a channel and show spikes
chIndx=19;
spikeTimes=spikes{chIndx};
spikeTrain=zeros(nSamples,1);
spikeTrain(spikeTimes)=1;
chWaveforms=waveform{chIndx};
figure;
subplot(211); hold on
plot(data(chIndx,:));
yl=ylim;
stem(1:nSamples,yl(2)*spikeTrain);
xlim([1 60*fs]); % show the first minute only
subplot(223)
plot(chWaveforms');
subplot(224);
plot(mean(chWaveforms,1));

% %%
% chidx=10;
% nfft=size(data,2);
% A=fft(data,nfft,2);
% freqs=(0:nfft-1)/nfft*fs;
% 
% figure
% subplot(221);
% plot(freqs,abs(A(chidx,:))); xlim([0 1000]);
% subplot(222)
% plot(data(chidx,:));
% 
% %%
% 
% 
% %%
% sampleStart=randi(size(data,2),1);
% sampleEnd=sampleStart+0.1*fs;
% dataShow=data(:,sampleStart:sampleEnd);
% time=(-20:19)/fs;
% figure
% for ch=1:32
% [sems,mus]=nansem(waveform{ch},1);
% subplot(4,8,ch);
% shadedErrorBar( time ,mus,sems);
% end
% figure
% for ch=1:32
% subplot(4,8,ch);
% plot(dataShow(ch,:));
% end
% 
% %%
% ch2show=15;
% [sems,mus]=nansem(waveform{ch},1);
% figure
% subplot(221); 
% plot((0:size(dataShow,2)-1)/fs,dataShow(ch,:),'k');
% axis tight
% ylabel('Amplitude (\muV)');
% xlabel('Time (s)');
% subplot(222);
% shadedErrorBar( time*1000 ,mus,sems);
% xlabel('Time (ms)');
% ylabel('Amplitude (\muV)')
% print -dpng ../figures/rpprMUA
% crop('../figures/rpprMUA.png',0);
