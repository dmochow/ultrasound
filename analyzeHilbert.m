clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));

%%
rStrs={'r012','r013','r014','r015','r016','r017','r018','r019'};
orders={ '40a','40s','10s','10a' ; ...
    '10a','10s','40a','40s'; ...
    '40s','40a','10s','10a'; ...
    '10s','10a','40a','40s'; ...
    '40s','40a','10a','10s'; ...
    '10a','10s','40s','40a';...
    '40a','40s','10a','10s';...
    '10s','10a','40s','40a'};

% excluding r018
% rStrs={'r012','r013','r014','r015','r016','r017','r019'};
% orders={ '40a','40s','10s','10a' ; ...
%     '10a','10s','40a','40s'; ...
%     '40s','40a','10s','10a'; ...
%     '10s','10a','40a','40s'; ...
%     '40s','40a','10a','10s'; ...
%     '10a','10s','40s','40a';...
%     '10s','10a','40s','40a'};



fs=500; % rate at which envelopes are sampled
fsd=50; % envelopes were downsampled to this
caChan=[21 12 22 11 23]; % hippocampus
icChans=[33:42]; %ICs
timeBsl=[-10*60*fsd:1:-1]/fsd/60;
timePost=[3*60*fsd+1:1:33*60*fsd]/fsd/60;
timeCat=[timeBsl timePost];
titles={'Active 40 Hz','Sham 40 Hz','Active 10 Hz','Sham 10 Hz'};
xl=[-2 13];

%%
for r=1:numel(rStrs)
    
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

grandMeans=squeeze(mean(allData,2));
grandMeansCorr=grandMeans;
% try to regress out extrapolated baseline
for c=1:8
    bslx=timeBsl(:);
    bsly=grandMeans(c,1:numel(timeBsl)); bsly=bsly(:); 
    B=regress(bsly,[bslx ones(numel(bslx),1)]);
    posty=B(2)+timePost*B(1);
    bslcorry=B(2)+timeBsl*B(1);
    grandMeansCorr(c,1:numel(timeBsl))=grandMeans(c,1:numel(timeBsl))-bslcorry;
    grandMeansCorr(c,numel(timeBsl)+1:end)=grandMeans(c,numel(timeBsl)+1:end)-posty;
end

% uncorrected
figure; 
for c=1:8
hs(c)=subplot(2,4,c); hold on
plot(timeBsl,grandMeans(c,1:numel(timeBsl)))
plot(timePost,grandMeans(c,numel(timeBsl)+1:end))
if c<5, htit(c)=title(titles{c});  end
xlim(xl);
end
set(get(hs(1),'Ylabel'),'String','40 Hz Power');
set(get(hs(5),'Ylabel'),'String','10 Hz Power');
set(hs(1:4),'Ylim',[0 70]);
set(hs(5:8),'Ylim',[0 200]);
%set(hs(1:8),'Xlim',[-10 30]);
set(hs(1:8),'Xlim',xl);
print -dpng ../figures/envelopes_allRats_

% baseline corrected
figure; 
for c=1:8
hs(c)=subplot(2,4,c); hold on
plot(timeBsl,grandMeansCorr(c,1:numel(timeBsl)))
plot(timePost,grandMeansCorr(c,numel(timeBsl)+1:end))
if c<5, htit(c)=title(titles{c});  end
end
set(get(hs(1),'Ylabel'),'String','40 Hz Power');
set(get(hs(5),'Ylabel'),'String','10 Hz Power');
set(hs(1:4),'Ylim',[0 70]);
set(hs(5:8),'Ylim',[0 200]);
%set(hs(1:8),'Xlim',[-10 30]);
set(hs(1:8),'Xlim',xl);
print -dpng ../figures/corrected_envelopes_allRats_


%%
% baseline corrected, active - sham
figure; 
for c=1:4
hs2(c)=subplot(2,2,c); hold on
plot(timeBsl,   grandMeansCorr(c*2-1,1:numel(timeBsl))  - grandMeansCorr(c*2,1:numel(timeBsl))     )
plot(timePost,  grandMeansCorr(c*2-1,numel(timeBsl)+1:end)  - grandMeansCorr(c*2,numel(timeBsl)+1:end))
end
set(get(hs2(1),'Title'),'String','Active 40 - Sham 40');
set(get(hs2(2),'Title'),'String','Active 10 - Sham 10');
set(get(hs2(1),'Ylabel'),'String','40 Hz Power');
set(get(hs2(3),'Ylabel'),'String','10 Hz Power');
set(hs2(1:2),'Ylim',[-40 60]);
set(hs2(3:4),'Ylim',[-150 150]);
set(hs2(1:4),'Xlim',[-10 40]);
print -dpng ../figures/corrected_activeminussham_envelopes_allRats_


%%
% % stats
% for c=1:8
%     x=timePost'; x=[x ones(numel(x),1)];
%     y=grandMeans(c,numel(timeBsl)+1:end)';
%     B=regress(y,x);
%     beta(c)=B(1); beta0(c)=B(2);
% end
% 
% figure;
% subplot(221);
% bar(1:4,beta(1:4));
% set(gca,'XtickLabel',{'Active 40','Sham 40','Active 10','Sham 10'});
% subplot(222);
% bar(1:4,beta(5:8));
% set(gca,'XtickLabel',{'Active 40','Sham 40','Active 10','Sham 10'});
% subplot(223);
% bar(1:4,beta0(1:4));
% set(gca,'XtickLabel',{'Active 40','Sham 40','Active 10','Sham 10'});
% subplot(224);
% bar(1:4,beta0(5:8));
% set(gca,'XtickLabel',{'Active 40','Sham 40','Active 10','Sham 10'});



