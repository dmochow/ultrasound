clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));

%% 
rStr='r012';
order={'40a','40s','10s','10a'};
fs=500; % rate at which envelopes are sampled
fsd=50; % envelopes were downsampled to this
%% 
% grab data for STIMULATION AT 40 HZ
freqStr='40HZ';
activeFolder=['../data/' rStr '/active' freqStr(1:2) ];
filenames=dir(fullfile(activeFolder,['*' rStr '_40HZ_ACTIVE_hilbert.mat']));
activeFilename=filenames(1).name;
shamFolder=['../data/' rStr '/sham' freqStr(1:2) ];
filenames=dir(fullfile(shamFolder,['*' rStr '_40HZ_SHAM_hilbert.mat']));
shamFilename=filenames(1).name;
varsActiveStim40=load(fullfile(activeFolder,activeFilename));
varsShamStim40=load(fullfile(shamFolder,shamFilename));
%% 
% grab data for STIMULATION AT 10 HZ
freqStr='10HZ';
activeFolder=['../data/' rStr '/active' freqStr(1:2) ];
filenames=dir(fullfile(activeFolder,['*' rStr '_10HZ_ACTIVE_hilbert.mat']));
activeFilename=filenames(1).name;
shamFolder=['../data/' rStr '/sham' freqStr(1:2) ];
filenames=dir(fullfile(shamFolder,['*' rStr '_10HZ_SHAM_hilbert.mat']));
shamFilename=filenames(1).name;
varsActiveStim10=load(fullfile(activeFolder,activeFilename));
varsShamStim10=load(fullfile(shamFolder,shamFilename));
%%
caChan=[21 12 22 11 23]; % hippocampus
icChans=[33:42]; %ICs
timeBsl=[-10*60*fsd:1:-1]/fsd/60;
timePost=[3*60*fsd+1:1:33*60*fsd]/fsd/60;

% cat data 2 show
figure;

for i=1:4

    cstr=order{i};
    
    switch cstr
        case '40a'
            dataCat=mean(varsActiveStim40.env40ds(caChan,:),1);
            subplot(2,4,i); hold on
            plot(timeCat,dataCat);
            ylabel('40 Hz envelope');
            xlabel('Time (min.)');
            title('Active 40 Hz TUS');
        case '40s'
        case '10a'
        case '10s'
            
    end
    
end


    
timeCat=[timeBsl timePost];



dataCat=mean(varsShamStim40.env40ds(caChan,:),1);
subplot(222); hold on
plot(timeCat,dataCat);
xlabel('Time (min.)');
title('Sham 40 Hz TUS');

dataCat=mean(varsActiveStim10.env10ds(caChan,:),1);
subplot(223); hold on
plot(timeCat,dataCat);
ylabel('10 Hz envelope');
xlabel('Time (min.)');
title('Active 10 Hz TUS');

dataCat=mean(varsShamStim10.env10ds(caChan,:),1);
subplot(224); hold on
plot(timeCat,dataCat);
xlabel('Time (min.)');
title('Sham 10 Hz TUS');


end




% %%
% % baseline to pre-TUS period
% [nChannels,nSamples,nBlocks]=size(varsActiveStim40.env40ds);
% bslInds= nSamples-round(30*fs)+1:nSamples; % baseline period (here 30 secs)
% 
% % active 40 env 40
% bsl=repmat(mean(varsActiveStim40.env40ds(:,bslInds,1),2),[1 nSamples nBlocks]);
% bslActive40Env40=(varsActiveStim40.env40ds-bsl)./bsl;
% tmp=bslActive40Env40(:,:,2:end); % omit baseline
% delsActive40Env40=tmp(:,:);
% delsActive40Env40=delsActive40Env40(probeMap,:);
% 
% % sham 40  env 40
% bsl=repmat(mean(varsShamStim40.env40ds(:,bslInds,1),2),[1 nSamples nBlocks]);
% bslSham40Env40=(varsShamStim40.env40ds-bsl)./bsl;
% tmp=bslSham40Env40(:,:,2:end); % omit baseline
% delsSham40Env40=tmp(:,:);
% delsSham40Env40=delsSham40Env40(probeMap,:);
% 
% % active 40 env 10
% bsl=repmat(mean(varsActiveStim40.env10ds(:,bslInds,1),2),[1 nSamples nBlocks]);
% bslActive40Env10=(varsActiveStim40.env10ds-bsl)./bsl;
% tmp=bslActive40Env10(:,:,2:end); % omit baseline
% delsActive40Env10=tmp(:,:);
% delsActive40Env10=delsActive40Env10(probeMap,:);
% 
% % sham 40 env 10 
% bsl=repmat(mean(varsShamStim40.env10ds(:,bslInds,1),2),[1 nSamples nBlocks]);
% bslSham40Env10=(varsShamStim40.env10ds-bsl)./bsl;
% tmp=bslSham40Env10(:,:,2:end); % omit baseline
% delsSham40Env10=tmp(:,:);
% delsSham40Env10=delsSham40Env10(probeMap,:);
% 
% delDelStim40Env40=delsActive40Env40-delsSham40Env40;
% delDelStim40Env10=delsActive40Env10-delsSham40Env10;
% 
% 
% 
% 
% % baseline to pre-TUS period
% [nChannels,nSamples,nBlocks]=size(varsActiveStim10.env40ds);
% bslInds= nSamples-round(30*fs)+1:nSamples; % baseline period (here 30 secs)
% 
% % active 10 env 40
% bsl=repmat(mean(varsActiveStim10.env40ds(:,bslInds,1),2),[1 nSamples nBlocks]);
% bslActive10Env40=(varsActiveStim10.env40ds-bsl)./bsl;
% tmp=bslActive10Env40(:,:,2:end); % omit baseline
% delsActive10Env40=tmp(:,:);
% delsActive10Env40=delsActive10Env40(probeMap,:);
% 
% % sham 10  env 40
% bsl=repmat(mean(varsShamStim10.env40ds(:,bslInds,1),2),[1 nSamples nBlocks]);
% bslSham10Env40=(varsShamStim10.env40ds-bsl)./bsl;
% tmp=bslSham10Env40(:,:,2:end); % omit baseline
% delsSham10Env40=tmp(:,:);
% delsSham10Env40=delsSham10Env40(probeMap,:);
% 
% % active 10 env 10
% bsl=repmat(mean(varsActiveStim10.env10ds(:,bslInds,1),2),[1 nSamples nBlocks]);
% bslActive10Env10=(varsActiveStim10.env10ds-bsl)./bsl;
% tmp=bslActive10Env10(:,:,2:end); % omit baseline
% delsActive10Env10=tmp(:,:);
% delsActive10Env10=delsActive10Env10(probeMap,:);
% 
% % sham 10 env 10 
% bsl=repmat(mean(varsShamStim10.env10ds(:,bslInds,1),2),[1 nSamples nBlocks]);
% bslSham10Env10=(varsShamStim10.env10ds-bsl)./bsl;
% tmp=bslSham10Env10(:,:,2:end); % omit baseline
% delsSham10Env10=tmp(:,:);
% delsSham10Env10=delsSham10Env10(probeMap,:);
% 
% delDelStim10Env40=delsActive10Env40-delsSham10Env40;
% delDelStim10Env10=delsActive10Env10-delsSham10Env10;
% 
% 
% %%
% figure
% chInd=[21 12 22 11 23];
% figure; 
% subplot(221); hold on
% plot(mean(delDelStim40Env40(chInd,:),1));
% subplot(222);
% plot(mean(delDelStim40Env10(chInd,:),1));
% subplot(223);
% plot(mean(delDelStim10Env40(chInd,:),1));
% subplot(224);
% plot(mean(delDelStim10Env10(chInd,:),1));
% 
% %%
% time=(0:size(delsActive40Env40,2)-1)/size(delsActive40Env40,2)*30;
% plotInds=1:1:90000;
% xl=[0 30];
% figure
% chInd=[21 12 22 11 23];
% %chInd=[33];
% figure; 
% subplot(221); hold on
% plot(time(plotInds),mean(delsActive40Env40(chInd,plotInds),1));
% plot(time(plotInds),mean(delsSham40Env40(chInd,plotInds),1));
% xlim(xl);
% subplot(222); hold on
% plot(time,mean(delsActive40Env10(chInd,:),1));
% plot(time,mean(delsSham40Env10(chInd,:),1));
% xlim(xl);
% subplot(223); hold on
% plot(time,mean(delsActive10Env40(chInd,:),1));
% plot(time,mean(delsSham10Env40(chInd,:),1));
% xlim(xl);
% subplot(224); hold on
% plot(time,mean(delsActive10Env10(chInd,:),1));
% plot(time,mean(delsSham10Env10(chInd,:),1));
% xlim(xl);


