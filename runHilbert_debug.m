clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
%%
rStr='r019';
freqStr='40HZ';
condStr='active';
subfolderStr='active40';
saveDataFileStr=[rStr '_' freqStr '_' condStr '_hilbert'];
targetFreq=str2num(freqStr(1:2));
if targetFreq==10
    offTargetFreq=40;
else
    offTargetFreq=10;
end
dataFolder=['../data/' rStr '/' subfolderStr '/' ];

%%
rerefChannel=NaN;
fs=30000;
dsr=60; %downsampling ratio
flHz=1; % lowest frequency that we're interested in
fhHz=250; % highest frequency we interested in
butterOrder=2;
probeMap=[17,16,18,15,19,14,20,13,...
    21,12,22,11,23,10,24,9,...
    25,8,26,7,27,6,28,5,...
    29,4,30,3,31,2,32,1]';
filenames=dir([dataFolder '*.mat']);
nFiles=numel(filenames);
for f=1:nFiles
    tFilename=upper(filenames(f).name);
    if contains(tFilename,'BASELINE') || contains(tFilename,'PRE') 
        bslFilename=tFilename;
    elseif contains(tFilename,'POSTTUS_1') || contains(tFilename,'POS_1') 
        post1Filename=tFilename;
    elseif contains(tFilename,'POSTTUS_2') || contains(tFilename,'POS_2')
        post2Filename=tFilename;
    elseif contains(tFilename,'POSTTUS_3') || contains(tFilename,'POS_3')
        post3Filename=tFilename;
    else
    end
end
matFilenames{1}=fullfile(dataFolder,bslFilename);
matFilenames{2}=fullfile(dataFolder,post1Filename);
matFilenames{3}=fullfile(dataFolder,post2Filename);
matFilenames{4}=fullfile(dataFolder,post3Filename);
samples={60*0*fs+1:60*10*fs;
    60*0*fs+1:60*10*fs;
    60*0*fs+1:60*10*fs;
    60*0*fs+1:60*10*fs}; % 10 min



%%
nFilenames=numel(matFilenames);
for f=1:nFilenames
    f
    load(matFilenames{f});
    data(:,:,f)=amplifier_data(:,samples{f}); clear amplifier_data;
end

%%
% rereference
if ~isnan(rerefChannel)
data = data - repmat (data(rerefChannel,:,:),[size(data,1) 1 1 ]);
end

%%
fl=flHz/(fs/2);
fh=fhHz/(fs/2);
[b,a]=butter(butterOrder,[fl fh]);
data=filter(b,a,data,[],2);


%%
% downsample data (assumes lowpass filtering already done)
dataDown=data(:,1:dsr:end,:);
fsr=fs/dsr;

%%
clear data

%%
nHarm=floor(fhHz/60);
for h=1:nHarm
    h
    [bn,an]=butter(2,[h*60-2 h*60+2]/fsr*2,'stop');
    dataDown=filter(bn,an,dataDown,[],2);
end

%% robust pca 
for b=1:size(dataDown,3)
    tmp=squeeze(dataDown(:,:,b)).'; % time x space
    [A_hat E_hat iter] = inexact_alm_rpca(tmp);
    dataDown(:,:,b)=A_hat.';
end

% %%
% % remove gross artifacts
% P=4;Q=4; zero=1; xtent=round(fsr/10);
% for f=1:size(dataDown,3)
%     X=dataDown(:,:,f);
%     Xout = nanBadSamples(X,P,Q,zero,xtent);
%     dataDown(:,:,f)=Xout;
% end

%%
% run ICA ?
X=dataDown(:,:); % combined 40 minute block
%[U,S,V]=svd(X',0); figure; plot(diag(S),'ko');
r=10; % estimated from above line after observing RPCA output
[Zica, Wica, Tica, muica] = kICA(X,r);
Rxx=cov(X');
Aica=Rxx*Wica*inv(Wica'*Rxx*Wica); figure; stem(Aica(probeMap,1));

%%
% comp=3;
% [zspec,fspec,tspec]=spectrogram(Zica(comp,:),5*fsr,1*fsr,[],fsr); 
% figure; subplot(221); imagesc([tspec(1) tspec(end)],[fspec(1) fspec(end)],abs(zspec));ylim([0 70]);
% subplot(222); stem(Aica(probeMap,comp));

%%
% add ICs to the data as auxiliary channels
Z3=reshape(Zica,[r size(dataDown,2) size(dataDown,3)]);
dataDown=cat(1,dataDown,Z3);


%%
% bandpass, hilbert, abs, save

% design filters
%[b40,a40]=butter(2,[38 42]/fsr*2);
[b40,a40]=butter(2,[36 44]/fsr*2);
if ~isstable(b40,a40), error('filter not stable bro'); end
%[b10,a10]=butter(2,[8 12]/fsr*2);
[b10,a10]=butter(2,[6 14]/fsr*2);
if ~isstable(b10,a10), error('filter not stable bro'); end

%%
% apply filters
dataDown40=filter(b40,a40,dataDown,[],2);
dataDown10=filter(b10,a10,dataDown,[],2);

for bl=1:size(dataDown,3)
    tmp=(hilbert(dataDown40(:,:,bl).')).';
    env40(:,:,bl)=abs(tmp);
    
    tmp=(hilbert(dataDown10(:,:,bl).')).';
    env10(:,:,bl)=abs(tmp);
end

%%
% subsample envelope to 100 Hz sampling
[bd,ad]=butter(2,50/(fsr*2),'low');
if ~isstable(bd,ad), error('filter not stable bro'); end
env40d=filter(bd,ad,env40,[],2);
env10d=filter(bd,ad,env10,[],2);
env40ds=env40d(:,1:10:end,:);
env10ds=env10d(:,1:10:end,:);

% %%
% nfft=2^nextpow2(size(env40,2));
% freqs=(0:nfft-1)/nfft*fsr;
% E40=fft(env40d,nfft,2);
% E10=fft(env10d,nfft,2);
% figure;
% subplot(211);
% plot(freqs,abs(E40(:,:,1))); xlim([0.5 50]);
% subplot(212);
% plot(freqs,abs(E10(:,:,1))); xlim([0.5 50]);
% 
% %%
% 
% %figure
% %%
% s1=env40(22,:,1);
% s2=env40(22,:,2);
% figure; hold on; plot(s1); plot(s2);



%%
% save for input to a second script
save(fullfile(dataFolder,[date '-' saveDataFileStr]),'fsr','Wica','Aica','Tica','Rxx',...
    'env40ds','env10ds','rerefChannel','butterOrder','flHz','fhHz','probeMap');
    
    

% 
% %%
% clear PXX
% window=2*fsr;
% noverlap=0*fsr;
% winlen=2*fsr;
% winshift=2*fsr;
% nWins=floor(size(dataDown,2)/winlen);
% NW=4; % time-bandwidth product for multitaper method
% nfft=2^nextpow2(winlen);
% for f=1:size(dataDown,3)
%     x=dataDown(:,:,f);
%     x=x.'; % time x space
%     for w=1:nWins
%         windx=(w-1)*winshift+1:(w-1)*winshift+winlen;
%         [pxx,freqs]=pmtm(x(windx,:),NW,nfft,fsr); 
%         PXX(:,:,w,f)=pxx;
%     end
% end
% 
% 
% %%
% % use median to reduce effect of outliers
% mPXX=squeeze(median(PXX,3));
% mpsc=(mPXX(:,:,3)-mPXX(:,:,1))./mPXX(:,:,1)*100;
% 
% %%
% chIds=probeMap;
% nChans=numel(chIds);
% nRows=4; nCols=ceil(nChans/nRows);
% 
% %delAllS= ( allS(:,:,:,3)-allS(:,:,:,1) ) ./ allS(:,:,:,1) ;
% %PSC= ( PXX(:,:,3)-PXX(:,:,1) ) ./ PXX(:,:,1);
% 
% figure;
% for ch=1:nChans
%     hs(ch)=subplot(nRows,nCols,ch); hold on
%     plot(freqs,mpsc(:,chIds(ch))); xlim([0 50]); yl=ylim;
%     %pcolor(t,fr,delAllS(:,:,chIds(ch))); shading interp;
%     %axis square
%     %axis tight
%     %caxis([-5 5]); ylim([0 50]); colormap jet;
%     %plot([t(1) t(end)],[targetFreq targetFreq],'--k','LineWidth',1);
%     plot([targetFreq targetFreq],[yl(1) yl(2)],'--k','LineWidth',1);
%     title(['Ch. ' num2str(chIds(ch))]);
% end
% 
% set(get(hs(1,1),'YLabel'),'String','Percent Change');
% set(get(hs(1,1),'XLabel'),'String','Freq (Hz)');
% 
% print('-dpng',['../figures/PSC_' rStr '_' freqStr '_' condStr ]);
% 
% %%
% % temporal evolution at target frequency
% [~,freqIndx]=min(abs(targetFreq-freqs));
% spspTargetPre=squeeze(PXX(freqIndx,probeMap,:,1));
% spspTargetPost=squeeze(PXX(freqIndx,probeMap,:,3));
% [~,offFreqIndx]=min(abs(offTargetFreq-freqs));
% spspOffTargetPre=squeeze(PXX(offFreqIndx,probeMap,:,1));
% spspOffTargetPost=squeeze(PXX(offFreqIndx,probeMap,:,3));
% figure
% subplot(221);
% imagesc([spspTargetPre spspTargetPost]);
% subplot(222);
% imagesc([spspOffTargetPre spspOffTargetPost]);
% subplot(223); hold on
% plot([spspTargetPre spspTargetPost].')
% plot(mean([spspTargetPre spspTargetPost]),'-k','LineWidth',2)
% subplot(224); hold on
% plot([spspOffTargetPre spspOffTargetPost].')
% plot(mean([spspOffTargetPre spspOffTargetPost]),'-k','LineWidth',2)
% 
% %%
% save(fullfile(dataFolder,saveDataFileStr),'PXX','freqs','fsr','rerefChannel','probeMap');
