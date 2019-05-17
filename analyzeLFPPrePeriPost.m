clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
%%
rStr='r010';
freqStr='40HZ';
condStr='SHAM';
saveDataFileStr=[rStr '_' freqStr '_' condStr '_PXX'];
targetFreq=str2num(freqStr(1:2));
if targetFreq==10
    offTargetFreq=40;
else
    offTargetFreq=10;
end
dataFolder=['../data/' rStr '/' ];
rerefChannel=NaN;
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
    tFilename=upper(filenames(f).name);
    if contains(tFilename,'PRE') && contains(tFilename,freqStr) && contains(tFilename,condStr)
        preFilename=tFilename;
    elseif contains(tFilename,'DUR') && contains(tFilename,freqStr) && contains(tFilename,condStr)
        periFilename=tFilename;
    elseif contains(tFilename,'POS') && contains(tFilename,freqStr) && contains(tFilename,condStr)
        postFilename=tFilename;
    else
    end
end
matFilenames{1}=fullfile(dataFolder,preFilename);
matFilenames{2}=fullfile(dataFolder,periFilename);
matFilenames{3}=fullfile(dataFolder,postFilename);
samples={60*0*fs+1:60*3*fs;
    60*0*fs+1:60*3*fs;
    60*0*fs+1:60*3*fs}; % 3 min

% % only for r007
% samples={60*7*fs+1:60*10*fs;
%     60*0*fs+1:60*3*fs;
%     60*0*fs+1:60*3*fs}; % 3 min


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
nHarm=floor(fhHz/60);
for h=1:nHarm
    h
    [bn,an]=butter(2,[h*60-2 h*60+2]/fs*2,'stop');
    data=filter(bn,an,data,[],2);
end


%%
% downsample data (assumes lowpass filtering already done)
dataDown=data(:,1:dsr:end,:);
fsr=fs/dsr;

%%
P=4;Q=4; zero=1; xtent=round(fsr/10);
for f=1:size(dataDown,3)
    X=dataDown(:,:,f);
    Xout = nanBadSamples(X,P,Q,zero,xtent);
    dataDown(:,:,f)=Xout;
end

%%
clear PXX
window=2*fsr;
noverlap=0*fsr;
winlen=2*fsr;
winshift=2*fsr;
nWins=floor(size(dataDown,2)/winlen);
NW=4; % time-bandwidth product for multitaper method
nfft=2^nextpow2(winlen);
for f=1:size(dataDown,3)
    x=dataDown(:,:,f);
    x=x.'; % time x space
    for w=1:nWins
        windx=(w-1)*winshift+1:(w-1)*winshift+winlen;
        [pxx,freqs]=pmtm(x(windx,:),NW,nfft,fsr); 
        PXX(:,:,w,f)=pxx;
    end
end


%%
% use median to reduce effect of outliers
mPXX=squeeze(median(PXX,3));
mpsc=(mPXX(:,:,3)-mPXX(:,:,1))./mPXX(:,:,1)*100;

%%
chIds=probeMap;
nChans=numel(chIds);
nRows=4; nCols=ceil(nChans/nRows);

%delAllS= ( allS(:,:,:,3)-allS(:,:,:,1) ) ./ allS(:,:,:,1) ;
%PSC= ( PXX(:,:,3)-PXX(:,:,1) ) ./ PXX(:,:,1);

figure;
for ch=1:nChans
    hs(ch)=subplot(nRows,nCols,ch); hold on
    plot(freqs,mpsc(:,chIds(ch))); xlim([0 50]); yl=ylim;
    %pcolor(t,fr,delAllS(:,:,chIds(ch))); shading interp;
    %axis square
    %axis tight
    %caxis([-5 5]); ylim([0 50]); colormap jet;
    %plot([t(1) t(end)],[targetFreq targetFreq],'--k','LineWidth',1);
    plot([targetFreq targetFreq],[yl(1) yl(2)],'--k','LineWidth',1);
    title(['Ch. ' num2str(chIds(ch))]);
end

set(get(hs(1,1),'YLabel'),'String','Percent Change');
set(get(hs(1,1),'XLabel'),'String','Freq (Hz)');

print('-dpng',['../figures/PSC_' rStr '_' freqStr '_' condStr ]);

%%
% temporal evolution at target frequency
[~,freqIndx]=min(abs(targetFreq-freqs));
spspTargetPre=squeeze(PXX(freqIndx,probeMap,:,1));
spspTargetPost=squeeze(PXX(freqIndx,probeMap,:,3));
[~,offFreqIndx]=min(abs(offTargetFreq-freqs));
spspOffTargetPre=squeeze(PXX(offFreqIndx,probeMap,:,1));
spspOffTargetPost=squeeze(PXX(offFreqIndx,probeMap,:,3));
figure
subplot(221);
imagesc([spspTargetPre spspTargetPost]);
subplot(222);
imagesc([spspOffTargetPre spspOffTargetPost]);
subplot(223); hold on
plot([spspTargetPre spspTargetPost].')
plot(mean([spspTargetPre spspTargetPost]),'-k','LineWidth',2)
subplot(224); hold on
plot([spspOffTargetPre spspOffTargetPost].')
plot(mean([spspOffTargetPre spspOffTargetPost]),'-k','LineWidth',2)

%%
save(fullfile(dataFolder,saveDataFileStr),'PXX','freqs','fsr','rerefChannel','probeMap');
