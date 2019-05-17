clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));

%
rStr='r015';
dateStr='07-May-2019-';
fs=500; % rate at which envelopes are sampled

%% STIMULATION AT 40 HZ
freqStr='40HZ';
activeFolder=['../data/' rStr '/active' freqStr(1:2) ];
%activeFilename='26-Feb-2019-r012_40HZ_ACTIVE_hilbert.mat';
activeFilename=[dateStr rStr '_40HZ_ACTIVE_hilbert.mat'];
shamFolder=['../data/' rStr '/sham' freqStr(1:2) ];
%shamFilename='06-Mar-2019-r012_40HZ_SHAM_hilbert.mat';
shamFilename=[dateStr rStr '_40HZ_SHAM_hilbert.mat'];

varsActiveStim40=load(fullfile(activeFolder,activeFilename));
varsShamStim40=load(fullfile(shamFolder,shamFilename));

% baseline to pre-TUS period
[nChannels,nSamples,nBlocks]=size(varsActiveStim40.env40);
bslInds= nSamples-round(30*fs)+1:nSamples; % baseline period (here 30 secs)

% active 40 env 40
bsl=repmat(mean(varsActiveStim40.env40(:,bslInds,1),2),[1 nSamples nBlocks]);
bslActive40Env40=varsActiveStim40.env40-bsl;
tmp=bslActive40Env40(:,:,2:end); % omit baseline
delsActive40Env40=tmp(:,:);
delsActive40Env40=delsActive40Env40(varsActiveStim40.probeMap,:);

% sham 40  env 40
bsl=repmat(mean(varsShamStim40.env40(:,bslInds,1),2),[1 nSamples nBlocks]);
bslSham40Env40=varsShamStim40.env40-bsl;
tmp=bslSham40Env40(:,:,2:end); % omit baseline
delsSham40Env40=tmp(:,:);
delsSham40Env40=delsSham40Env40(varsShamStim40.probeMap,:);

% active 40 env 10
bsl=repmat(mean(varsActiveStim40.env10(:,bslInds,1),2),[1 nSamples nBlocks]);
bslActive40Env10=varsActiveStim40.env10-bsl;
tmp=bslActive40Env10(:,:,2:end); % omit baseline
delsActive40Env10=tmp(:,:);
delsActive40Env10=delsActive40Env10(varsActiveStim40.probeMap,:);

% sham 40 env 10 
bsl=repmat(mean(varsShamStim40.env10(:,bslInds,1),2),[1 nSamples nBlocks]);
bslSham40Env10=varsShamStim40.env10-bsl;
tmp=bslSham40Env10(:,:,2:end); % omit baseline
delsSham40Env10=tmp(:,:);
delsSham40Env10=delsSham40Env10(varsShamStim40.probeMap,:);

delDelStim40Env40=delsActive40Env40-delsSham40Env40;
delDelStim40Env10=delsActive40Env10-delsSham40Env10;


%% STIMULATION AT 10 HZ
freqStr='10HZ';
activeFolder=['../data/' rStr '/active' freqStr(1:2) ];
%activeFilename='07-May-2019-r012_10HZ_ACTIVE_hilbert.mat';
activeFilename=[dateStr rStr '_10HZ_ACTIVE_hilbert.mat'];
shamFolder=['../data/' rStr '/sham' freqStr(1:2) ];
%shamFilename='07-May-2019-r012_10HZ_SHAM_hilbert.mat';
shamFilename=[dateStr rStr '_10HZ_SHAM_hilbert.mat'];

varsActiveStim10=load(fullfile(activeFolder,activeFilename));
varsShamStim10=load(fullfile(shamFolder,shamFilename));

% baseline to pre-TUS period
[nChannels,nSamples,nBlocks]=size(varsActiveStim10.env40);
bslInds= nSamples-round(30*fs)+1:nSamples; % baseline period (here 30 secs)

% active 10 env 40
bsl=repmat(mean(varsActiveStim10.env40(:,bslInds,1),2),[1 nSamples nBlocks]);
bslActive10Env40=varsActiveStim10.env40-bsl;
tmp=bslActive10Env40(:,:,2:end); % omit baseline
delsActive10Env40=tmp(:,:);
delsActive10Env40=delsActive10Env40(varsActiveStim10.probeMap,:);

% sham 10  env 40
bsl=repmat(mean(varsShamStim10.env40(:,bslInds,1),2),[1 nSamples nBlocks]);
bslSham10Env40=varsShamStim10.env40-bsl;
tmp=bslSham10Env40(:,:,2:end); % omit baseline
delsSham10Env40=tmp(:,:);
delsSham10Env40=delsSham10Env40(varsShamStim10.probeMap,:);

% active 10 env 10
bsl=repmat(mean(varsActiveStim10.env10(:,bslInds,1),2),[1 nSamples nBlocks]);
bslActive10Env10=varsActiveStim10.env10-bsl;
tmp=bslActive10Env10(:,:,2:end); % omit baseline
delsActive10Env10=tmp(:,:);
delsActive10Env10=delsActive10Env10(varsActiveStim40.probeMap,:);

% sham 10 env 10 
bsl=repmat(mean(varsShamStim10.env10(:,bslInds,1),2),[1 nSamples nBlocks]);
bslSham10Env10=varsShamStim10.env10-bsl;
tmp=bslSham10Env10(:,:,2:end); % omit baseline
delsSham10Env10=tmp(:,:);
delsSham10Env10=delsSham10Env10(varsShamStim10.probeMap,:);

delDelStim10Env40=delsActive10Env40-delsSham10Env40;
delDelStim10Env10=delsActive10Env10-delsSham10Env10;


%%
figure
chInd=[21 12 22 11 23];
figure; 
subplot(221);
plot(mean(delDelStim40Env40(chInd,:),1));
subplot(222);
plot(mean(delDelStim40Env10(chInd,:),1));
subplot(223);
plot(mean(delDelStim10Env40(chInd,:),1));
subplot(224);
plot(mean(delDelStim10Env10(chInd,:),1));


% %tmp=varsActive.env40(:,:,2:4)-varsActive.env40(:,:,1);
% 
% 
% %dataFolder=['../data/' rStr '/' subfolderStr '/' ];
% %load(fullfile(dataFolder,[dateStr '-' saveDataFileStr]));
% %figFileStr=[date '-' rStr '-' freqStr '-' condStr '-envelope-time-series.png' ]
% 
% %%
% % downsample envelopes to 20 Hz sampling rate
% fsfinal=20;
% dsr=fsr/fsfinal;
% env40down=env40(:,1:dsr:end,:);
% env10down=env10(:,1:dsr:end,:);
% 
% %%
% env40down2d=env40down(:,:);
% env10down2d=env10down(:,:);
% 
% %%
% nChannels=size(env40down2d,1);
% nSamples=size(env40down2d,2);
% nRows=nChannels;
% nCols=2;
% time=(0:nSamples-1)/fsfinal;
% 
% figure
% for ch=1:nChannels
%     subplot(nRows,nCols,(ch-1)*2+1); hold on
%     plot(time,env40down2d(probeMap(ch),:));
%     yl=ylim;
%     plot([10*60 10*60],[yl(1) yl(2)],'--k');
%     
%     subplot(nRows,nCols,ch*2); hold on
%     plot(time,env10down2d(probeMap(ch),:));
%     yl=ylim;
%     plot([10*60 10*60],[yl(1) yl(2)],'--k');
% end
% 
% %%
% print('-dpng',fullfile('../figures/',figFileStr));
% crop(fullfile('../figures/',figFileStr),0);
% %%
% % % display by blocks
% % nBlocks=size(env40down,3);
% % figure
% % for bl=1:nBlocks
% %     subplot(4,1,bl); hold on
% %     plot(mean(env40down(:,:,bl),1))
% %     plot(mean(env10down(:,:,bl),1))
% % end