clear all; close all; clc
load ../data/smartboxFirstData; fs=30000;
%load ../data/iso3fixedground; fs=20000;


data=amplifier_data;
channelsKeep=1:32; % subset of channels that we want to keep
flHz=300; % lowest frequency that we're interested in
fhHz=3000; % highest frequency
fsr=floor(fs/fhHz);
fsd=fs/fsr;
butterOrder=4;
dimFilter=2;
fl=flHz/(fs/2);
fh=fhHz/(fs/2);
[b,a]=butter(butterOrder,[fl fh]);
data=filter(b,a,data,[],2);

nHarm=50;
for h=1:nHarm
    [bn,an]=butter(2,[h*60-2 h*60+2]/fs*2,'stop');
    data=filter(bn,an,data,[],2);
end

%data=downsample(data.',fsr).';

%% spike detection parameters
% spike detection threshold (multiple of median(abs)/0.6745)
threshStds=5;









%%
chidx=10;
nfft=size(data,2);
A=fft(data,nfft,2);
freqs=(0:nfft-1)/nfft*fs;

figure
subplot(221);
plot(freqs,abs(A(chidx,:))); xlim([0 1000]);
subplot(222)
plot(data(chidx,:));

%%
[spikes,waveform,mua,thresh] = detectSpikes(data(:,:),5,fs,'neg');

%%
sampleStart=randi(size(data,2),1);
sampleEnd=sampleStart+0.1*fs;
dataShow=data(:,sampleStart:sampleEnd);
time=(-20:19)/fs;
figure
for ch=1:32
[sems,mus]=nansem(waveform{ch},1);
subplot(4,8,ch);
shadedErrorBar( time ,mus,sems);
end
figure
for ch=1:32
subplot(4,8,ch);
plot(dataShow(ch,:));
end

%%
ch2show=15;
[sems,mus]=nansem(waveform{ch},1);
figure
subplot(221); 
plot((0:size(dataShow,2)-1)/fs,dataShow(ch,:),'k');
axis tight
ylabel('Amplitude (\muV)');
xlabel('Time (s)');
subplot(222);
shadedErrorBar( time*1000 ,mus,sems);
xlabel('Time (ms)');
ylabel('Amplitude (\muV)')
print -dpng ../figures/rpprMUA
crop('../figures/rpprMUA.png',0);
