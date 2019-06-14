% 06.05.19
% analyze envelopes of rats for which we have data collected during
% stimulation, where the probe was not moved in and out as before
% (analyzeHilbert.m)
clear all; close all; clc


%%
rStrs={'r022'};
orders={'40s','40a','10s','10a'};

fs=500; % rate at which envelopes are sampled
fsd=50; % envelopes were downsampled to this
%caChan=[21 12 22 11 23]; % hippocampus
caChan=[9 25 8 26]; % hippocampus
timeBsl=[-10*60*fsd:1:-1]/fsd/60;
timeStim=[0*fsd+1:1:3*60*fsd]/fsd/60;
timePost=[3*60*fsd+1:1:33*60*fsd]/fsd/60;
timeCat=[timeBsl timeStim timePost];
titles={'Active 40 Hz','Sham 40 Hz','Active 10 Hz','Sham 10 Hz'};
xl=[-2 13];
nRats=numel(rStrs);

%%
for r=1:nRats
    
    rStr=rStrs{r};
    %order=orders{r,:}

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
    
    allData(1,r,:)=mean(varsActiveStim40.env40ds(caChan,:),1);
    allData(2,r,:)=mean(varsShamStim40.env40ds(caChan,:),1);
    allData(3,r,:)=mean(varsActiveStim10.env40ds(caChan,:),1);
    allData(4,r,:)=mean(varsShamStim10.env40ds(caChan,:),1);
        
    allData(5,r,:)=mean(varsActiveStim40.env10ds(caChan,:),1);    
    allData(6,r,:)=mean(varsShamStim40.env10ds(caChan,:),1);
    allData(7,r,:)=mean(varsActiveStim10.env10ds(caChan,:),1);
    allData(8,r,:)=mean(varsShamStim10.env10ds(caChan,:),1);
    
    
end

%%
% regress out extrapolated baseline from individual recordings
allDataCorr=allData;
bslx=timeBsl(:);
for c=1:8
    for r=1:nRats
        bsly=squeeze(allData(c,r,1:numel(timeBsl)));
        bsly=bsly(:);
        B=regress(bsly,[bslx ones(numel(bslx),1)]);
        
        posty=B(2)+timePost*B(1);
        stimy=B(2)+timeStim*B(1);
        bslcorry=B(2)+timeBsl*B(1);
        
        %
        inp=squeeze(allData(c,r,1:numel(timeBsl))).';
        otp=inp-bslcorry;
        allDataCorr(c,r,1:numel(timeBsl))=otp;
        
        %
        inp=squeeze(allData(c,r,numel(timeBsl)+1:numel(timeBsl)+numel(timeStim))).';
        otp=inp-stimy;
        allDataCorr(c,r,numel(timeBsl)+1:numel(timeBsl)+numel(timeStim))=otp;
        
        %
        inp=squeeze(allData(c,r,numel(timeBsl)+numel(timeStim)+1:end)).';
        otp=inp-posty;
        allDataCorr(c,r,numel(timeBsl)+numel(timeStim)+1:end)=otp;
    end
end
% allDataCorr=allData; % omit regressing of the baseline trend



%%
dsr=50;
grandMeansCorr=squeeze(mean(allDataCorr,2));
grandSemsCorr=squeeze(std(allDataCorr,[],2))/sqrt(nRats);
for c=1:8
    c
    grandMeansCorrDown(c,:)=resample(grandMeansCorr(c,:),1,dsr);
    grandSemsCorrDown(c,:)=resample(grandSemsCorr(c,:),1,dsr);
end
timeBslDown=timeBsl(1:dsr:end); 
timeStimDown=timeStim(1:dsr:end); 
timePostDown=timePost(1:dsr:end); 
timeDown=[timeBslDown timeStimDown timePostDown];

% grandMeansCorrDown=grandMeansCorr;
% timeBslDown=timeBsl;
% timePostDown=timePost;

%%
% show subject averaged downsampled curve
conds2show=[1 2 5 6];
figure; 
hs(1)=subplot(221); hold on
hsh(1)=plot(timeDown, grandMeansCorrDown(1,:));
hsh(2)=plot(timeDown, grandMeansCorrDown(2,:));
ylabel('40 Hz power');
hlg=legend('active 40', 'sham 40');
%xlim([0 3]);

hs(2)=subplot(222); hold on
hsh(3)=plot(timeDown, grandMeansCorrDown(3,:));
hsh(4)=plot(timeDown, grandMeansCorrDown(4,:));
ylabel('40 Hz power');
hlg=legend('active 10', 'sham 10');
%xlim([0 3]);

hs(3)=subplot(223); hold on
hsh(5)=plot(timeDown, grandMeansCorrDown(5,:));
hsh(6)=plot(timeDown, grandMeansCorrDown(6,:));
ylabel('10 Hz power');
hlg=legend('active 40', 'sham 40');
%xlim([0 3]);

hs(4)=subplot(224); hold on
hsh(7)=plot(timeDown, grandMeansCorrDown(7,:));
hsh(8)=plot(timeDown, grandMeansCorrDown(8,:));
ylabel('10 Hz power');
hlg=legend('active 10', 'sham 10');
%xlim([0 3]);







%hs(2)=subplot(212); hold on
%plot(timeBslDown,grandMeansCorrDown(5,1:numel(timeBslDown)))
%plot(timePostDown,grandMeansCorrDown(5,numel(timeBslDown)+1:end))
%plot(timeBslDown,grandMeansCorrDown(6,1:numel(timeBslDown)),'--')
%plot(timePostDown,grandMeansCorrDown(6,numel(timeBslDown)+1:end),'--')
%set(get(hs(2),'Ylabel'),'String','10 Hz Power');

%set(hs(1:4),'Ylim',[0 70]);
%set(hs(5:8),'Ylim',[0 200]);
%set(hs(1:8),'Xlim',[-10 30]);
%set(hs(1:8),'Xlim',xl);
% print -dpng ../figures/mean_envelopes_allRats_








