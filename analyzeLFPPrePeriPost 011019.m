clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
%%
rStr='r010';
freqStr='40Hz';
targetFreq=10;
condStr='Sham';
dataFolder=['../data/' rStr '/' ];
rerefChannel=1;
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
samples={60*0*fs+1:60*3*fs;
    60*0*fs+1:60*3*fs;
    60*0*fs+1:60*3*fs}; % 3 min


%%
nFilenames=numel(matFilenames);
for f=1:nFilenames
    f
    load(matFilenames{f});
    data(:,:,f)=amplifier_data(:,samples{f}); clear amplifier_data;
end

%%
% rereference
data = data - repmat (data(rerefChannel,:,:),[size(data,1) 1 1 ]);

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
for f=1:size(dataDown,3)
    for ch=1:size(dataDown,1)
        [f,ch]
        x=dataDown(ch,:,f);
        [S,fr,t] = spectrogram(x,window,noverlap,[],fsr);
        allS(:,:,ch,f)=abs(S);
    end
    %     x=x.'; % timexspace
    %     [pxx,fr] = pwelch(x,window,noverlap,[],fsr);
    %     %pxx=pmtm(x); %fr=(0:size(PXX,1)-1)/size(x,1)*fsr;
    %     PXX(:,:,f)=pxx;
end


%%
chIds=probeMap;
nChans=numel(chIds);
nRows=4; nCols=ceil(nChans/nRows);

delAllS= ( allS(:,:,:,3)-allS(:,:,:,1) ) ./ allS(:,:,:,1) ;
%PSC= ( PXX(:,:,3)-PXX(:,:,1) ) ./ PXX(:,:,1);

figure;
for ch=1:nChans
    hs(ch)=subplot(nRows,nCols,ch); hold on
    %plot(fr,PSC(:,chIds(ch))); xlim([0 50]);
    pcolor(t,fr,delAllS(:,:,chIds(ch))); shading interp;
    axis square
    axis tight
    caxis([-5 5]); ylim([0 50]); colormap jet;
    plot([t(1) t(end)],[targetFreq targetFreq],'--k','LineWidth',1);
    title(['Ch. ' num2str(chIds(ch))]);
end

%set(get(hs(1,1),'YLabel'),'String','Percent Change');
%set(get(hs(1,1),'XLabel'),'String','Freq (Hz)');

print('-dpng',['../figures/delSpectrogram_' rStr '_' freqStr '_' condStr ]);

%%
% use median to reduce effect of outliers
meds=squeeze(median(allS,2));
pscMed=(meds(:,:,3)-meds(:,:,1))./meds(:,:,1);
figure;
plot(fr,pscMed(:,8));
