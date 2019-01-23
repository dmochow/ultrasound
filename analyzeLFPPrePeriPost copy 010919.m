clear all; close all; clc

%%
rStr='r010';
freqStr='10Hz';
condStr='Active';
%dataFolder=['../data/' rStr '/*.mat' ];
dataFolder=['../data/' rStr '/' ];
fs=30000;
dsr=30; %downsampling ratio
flHz=1; % lowest frequency that we're interested in
fhHz=500; % highest frequency we interested in
butterOrder=3;
probeMap=[17,16,18,15,19,14,20,13,...
    21,12,22,11,23,10,24,9,...
    25,8,26,7,27,6,28,5,...
    29,4,30,3,31,2,32,1]';
filenames=dir([dataFolder '*.mat']);
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
% rereference

data = data - repmat (data(32,:,:),[size(data,1) 1 1 ]);

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

%%
% downsample data (assumes lpf already done)
dataDown=data(:,1:dsr:end,:);
fsr=fs/dsr;

%%
clear PXX
window=2*fsr;
noverlap=1*fsr;
for f=1:3
    f
    x=dataDown(:,:,f);
    x=x.'; % timexspace
    [pxx,fr] = pwelch(x,window,noverlap,[],fsr);
    %[pxx,fr] = periodogram(x,[],[],fsr);
    PXX(:,:,f)=pxx;
end

%nfft=2^nextpow2(size(dataDown,2));
%D=fft(dataDown,nfft,2);
%freqs=(0:nfft-1)/nfft*fsr;

%%
% D1=D(:,:,1);
% D2=D(:,:,2);
% D3=D(:,:,3);
figure;
subplot(311);
plot(fr,mean(PXX(:,:,1),2)); xlim([0 200]);
subplot(312);
plot(fr,mean(PXX(:,:,2),2)); xlim([0 200]);
subplot(313);
plot(fr,mean(PXX(:,:,3),2)); xlim([0 200]);

%%
chIds=[26 25 8 1 32];
nChans=numel(chIds);
nRows=3; nCols=nChans;

figure;
for ch=1:nChans
    
    hs(1,ch)=subplot(nRows,nCols,ch);
    plot(fr,PXX(:,chIds(ch),1)); xlim([0 50]); title(['Pre Chan.' num2str(ch)])
    
    hs(2,ch)=subplot(nRows,nCols,ch+nCols);
    plot(fr,PXX(:,chIds(ch),2)); xlim([0 50]); title(['TUS Chan.' num2str(ch)])
    
    hs(3,ch)=subplot(nRows,nCols,ch+2*nCols);
    plot(fr,PXX(:,chIds(ch),3)); xlim([0 50]); title(['Post Chan.' num2str(ch)])
end

% set(hs(1,:),'Ylim',[0 600]);
% set(hs(2,:),'Ylim',[0 60000]);
% set(hs(3,:),'Ylim',[0 600]);
set(get(hs(1,1),'YLabel'),'String','Power Spectrum');
set(get(hs(1,1),'XLabel'),'String','Freq (Hz)');

% print -dpng ../figures/r009_LFP_prePeriPost
% crop('../figures/r009_LFP_prePeriPost.png',0);


%%
f10=22;
lfpwr=squeeze(mean(  PXX(1:f10-1,:,:),1  ));
lfpwratio=lfpwr(:,3)./lfpwr(:,1);

pwr10=squeeze(PXX(f10,:,:));
pwr10ratio=pwr10(:,3)./pwr10(:,1);
figure; stem(pwr10ratio);

% tuspwr=squeeze(PXX(f10,:,2));
% figure; scatter(tuspwr,lfpwratio);
% [sortvals,sortind]=sort(lfpwratio,'ascend');
% [sortvals2,sortind2]=sort(tuspwr,'descend');

%%

%[spikes,waveform,mua,thresh] = filterDetect(amplifier_data);

%%
% spectrograms
chIdx=1;
chDataPre=dataDown(chIdx,:,1);
chDataPeri=dataDown(chIdx,:,2);
chDataPost=dataDown(chIdx,:,3);
windows=fsr;
noverlap=fsr/2'
[S1,f,t]=spectrogram(chDataPre,window,noverlap,[],fsr);
[S2,f,t]=spectrogram(chDataPeri,window,noverlap,[],fsr);
[S3,f,t]=spectrogram(chDataPost,window,noverlap,[],fsr);
figure;
subplot(131);
pcolor(t,f,abs(S1)); shading interp; ylim([0 50]);
cax=caxis;
subplot(132);
pcolor(t,f,abs(S2)); shading interp;  ylim([0 50]);caxis(cax);
subplot(133);
pcolor(t,f,abs(S3)); shading interp;  ylim([0 50]);caxis(cax);