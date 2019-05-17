% LFP group analysis
clear all; close all; clc
rStr='r010';
freqStr='40Hz';
targetFreq=str2num(freqStr(1:2));
activeDataFileStr=[rStr '_' freqStr '_active_PXX'];
shamDataFileStr=[rStr '_' freqStr '_sham_PXX'];
dataFolder=['../data/' rStr '/' ];

%%
if targetFreq==10
    offTargetFreq=40;
else
    offTargetFreq=10;
end
load(fullfile(dataFolder,activeDataFileStr),'PXX','freqs','fsr','rerefChannel','probeMap');
aPXX=PXX;
load(fullfile(dataFolder,shamDataFileStr),'PXX','freqs','fsr','rerefChannel','probeMap');
sPXX=PXX;

%%
[~,targetFreqIndx]=min(abs(targetFreq-freqs));
[~,offTargetFreqIndx]=min(abs(offTargetFreq-freqs));

tsActiveTargetPre=squeeze(aPXX(targetFreqIndx,probeMap,:,1));
tsActiveTargetPost=squeeze(aPXX(targetFreqIndx,probeMap,:,3));

tsActiveOffTargetPre=squeeze(aPXX(offTargetFreqIndx,probeMap,:,1));
tsActiveOffTargetPost=squeeze(aPXX(offTargetFreqIndx,probeMap,:,3));

tsShamTargetPre=squeeze(sPXX(targetFreqIndx,probeMap,:,1));
tsShamTargetPost=squeeze(sPXX(targetFreqIndx,probeMap,:,3));

tsShamOffTargetPre=squeeze(sPXX(offTargetFreqIndx,probeMap,:,1));
tsShamOffTargetPost=squeeze(sPXX(offTargetFreqIndx,probeMap,:,3));

%%
figure
subplot(221); hold on
plot([tsActiveTargetPre tsActiveTargetPost].');
plot(mean([tsActiveTargetPre tsActiveTargetPost]),'-k','LineWidth',2)
title('Active On-Target');
subplot(222); hold on
plot([tsActiveOffTargetPre tsActiveOffTargetPost].');
plot(mean([tsActiveOffTargetPre tsActiveOffTargetPost]),'-k','LineWidth',2)
title('Active Off-Target');
subplot(223); hold on
plot([tsShamTargetPre tsShamTargetPost].');
plot(mean([tsShamTargetPre tsShamTargetPost]),'-k','LineWidth',2)
title('Sham On-Target');
subplot(224); hold on
plot([tsShamOffTargetPre tsShamOffTargetPost].');
plot(mean([tsShamOffTargetPre tsShamOffTargetPost]),'-k','LineWidth',2)
title('Sham Off-Target');



%%
amPXX=squeeze(median(aPXX,3));
ampsc=(amPXX(:,:,3)-amPXX(:,:,1))./amPXX(:,:,1)*100;
smPXX=squeeze(median(sPXX,3));
smpsc=(smPXX(:,:,3)-smPXX(:,:,1))./smPXX(:,:,1)*100;

% %%
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
%     plot(freqs,ampsc(:,chIds(ch)) - smpsc(:,chIds(ch)) ); xlim([0 50]); yl=ylim;
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