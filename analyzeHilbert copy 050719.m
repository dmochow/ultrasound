clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));

%%
rStr='r012';
freqStr='40HZ';
condStr='ACTIVE';
subfolderStr='sham40';
%dateStr='26-Feb-2019';
dateStr='06-Mar-2019';
saveDataFileStr=[rStr '_' freqStr '_' condStr '_hilbert'];
dataFolder=['../data/' rStr '/' subfolderStr '/' ];
load(fullfile(dataFolder,[dateStr '-' saveDataFileStr]));
figFileStr=[date '-' rStr '-' freqStr '-' condStr '-envelope-time-series.png' ]

%%
% downsample envelopes to 20 Hz sampling rate
fsfinal=20;
dsr=fsr/fsfinal;
env40down=env40(:,1:dsr:end,:);
env10down=env10(:,1:dsr:end,:);

%%
env40down2d=env40down(:,:);
env10down2d=env10down(:,:);

%%
nChannels=size(env40down2d,1);
nSamples=size(env40down2d,2);
nRows=nChannels;
nCols=2;
time=(0:nSamples-1)/fsfinal;

figure
for ch=1:nChannels
    subplot(nRows,nCols,(ch-1)*2+1); hold on
    plot(time,env40down2d(probeMap(ch),:));
    yl=ylim;
    plot([10*60 10*60],[yl(1) yl(2)],'--k');
    
    subplot(nRows,nCols,ch*2); hold on
    plot(time,env10down2d(probeMap(ch),:));
    yl=ylim;
    plot([10*60 10*60],[yl(1) yl(2)],'--k');
end

%%
print('-dpng',fullfile('../figures/',figFileStr));
crop(fullfile('../figures/',figFileStr),0);
%%
% % display by blocks
% nBlocks=size(env40down,3);
% figure
% for bl=1:nBlocks
%     subplot(4,1,bl); hold on
%     plot(mean(env40down(:,:,bl),1))
%     plot(mean(env10down(:,:,bl),1))
% end