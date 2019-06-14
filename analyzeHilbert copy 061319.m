% 05.24.19
% here we are correcting individual rats, then averaging
% in previous analyzeHibert, we were averaging and then correcting

clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));

%%
% rStrs={'r012','r013','r014','r015','r016','r017','r018','r019','r022'};
% orders={ '40a','40s','10s','10a' ; ...
%     '10a','10s','40a','40s'; ...
%     '40s','40a','10s','10a'; ...
%     '10s','10a','40a','40s'; ...
%     '40s','40a','10a','10s'; ...
%     '10a','10s','40s','40a';...
%     '40a','40s','10a','10s';...
%     '10s','10a','40s','40a';...
%     '40s','40a','10s','10a'};

%excluding r012
% rStrs={'r013','r014','r015','r016','r017','r018','r019','r022'};
% orders={'10a','10s','40a','40s'; ...
%     '40s','40a','10s','10a'; ...
%     '10s','10a','40a','40s'; ...
%     '40s','40a','10a','10s'; ...
%     '10a','10s','40s','40a';...
%     '40a','40s','10a','10s';...
%     '10s','10a','40s','40a';...
%      '40s','40a','10s','10a'};

%excluding r018
rStrs={'r012','r013','r014','r015','r016','r017','r019','r022','r023'};
orders={ '40a','40s','10s','10a' ; ...
    '10a','10s','40a','40s'; ...
    '40s','40a','10s','10a'; ...
    '10s','10a','40a','40s'; ...
    '40s','40a','10a','10s'; ...
    '10a','10s','40s','40a';...
    '10s','10a','40s','40a';...
    '40s','40a','10s','10a';...
    '10s','10a','40s','40a'};

%excluding r018
rStrs={'r012','r013','r014','r015','r016','r017','r019'};
orders={ '40a','40s','10s','10a' ; ...
    '10a','10s','40a','40s'; ...
    '40s','40a','10s','10a'; ...
    '10s','10a','40a','40s'; ...
    '40s','40a','10a','10s'; ...
    '10a','10s','40s','40a';...
    '10s','10a','40s','40a'};

% % angled
% rStrs={'r022','r023'};
% orders={ '40s','40a','10s','10a';...
%     '10s','10a','40s','40a'};
% 
% rStrs={'r022'};
% orders={ '40s','40a','10s','10a'};
% 
% rStrs={'r023'};
% orders={ '10s','10a','40s','40a'};





fs=500; % rate at which envelopes are sampled
fsd=50; % envelopes were downsampled to this
caChan1=[21 12 22 11 23]; % hippocampus for r012-r019
caChan2=[9 25 8 26]; % hippocampus for r022-
icChans=[33:42]; %ICs
timeBsl=[-10*60*fsd:1:-1]/fsd/60;
timePost=[3*60*fsd+1:1:33*60*fsd]/fsd/60;
timeCat=[timeBsl timePost];
titles={'Active 40 Hz','Sham 40 Hz','Active 10 Hz','Sham 10 Hz'};
xl=[-2 13];
nRats=numel(rStrs);

%%
for r=1:nRats
    
    rStr=rStrs{r};
    
    if str2double(rStr(2:end))>21
        caChan=caChan2;
    else
        caChan=caChan1;
    end
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
        bslcorry=B(2)+timeBsl*B(1);
        %
        inp=squeeze(allData(c,r,1:numel(timeBsl))).';
        otp=inp-bslcorry;
        %otp=inp;
        allDataCorr(c,r,1:numel(timeBsl))=otp;
        %
        inp=squeeze(allData(c,r,numel(timeBsl)+1:end)).';
        otp=inp-posty;
        %otp=inp;
        allDataCorr(c,r,numel(timeBsl)+1:end)=otp;
    end
end

% normalize
allDataCorr_z = allDataCorr./repmat(mean(allDataCorr.^2,3),[1 1 size(allDataCorr,3)]);

%%
dsr=50;
% un-normalized
grandMeansCorr=squeeze(mean(allDataCorr,2));
grandSemsCorr=squeeze(std(allDataCorr,[],2))/sqrt(nRats);

% normalized
% grandMeansCorr=squeeze(mean(allDataCorr_z,2));
% grandSemsCorr=squeeze(std(allDataCorr_z,[],2))/sqrt(nRats);
for c=1:8
    c
    grandMeansCorrDown(c,:)=resample(grandMeansCorr(c,:),1,dsr);
    grandSemsCorrDown(c,:)=resample(grandSemsCorr(c,:),1,dsr);
end
timeBslDown=timeBsl(1:dsr:end); 
timePostDown=timePost(1:dsr:end); 

% grandMeansCorrDown=grandMeansCorr;
% timeBslDown=timeBsl;
% timePostDown=timePost;

%%
% show subject averaged downsampled curve
figure; 
hs(1)=subplot(221); hold on
hsh(1)=shadedErrorBar(timeBslDown,grandMeansCorrDown(1,1:numel(timeBslDown))  , grandSemsCorrDown(1,1:numel(timeBslDown)) ,'b' ,1);
hsh(2)=shadedErrorBar(timePostDown,grandMeansCorrDown(1,numel(timeBslDown)+1:end),grandSemsCorrDown(1,numel(timeBslDown)+1:end) ,'b' ,1);
hsh(3)=shadedErrorBar(timeBslDown,grandMeansCorrDown(2,1:numel(timeBslDown))  , grandSemsCorrDown(2,1:numel(timeBslDown)) ,'r',1);
hsh(4)=shadedErrorBar(timePostDown,grandMeansCorrDown(2,numel(timeBslDown)+1:end),grandSemsCorrDown(2,numel(timeBslDown)+1:end),'r',1);
ylabel('40 Hz Power');

hs(2)=subplot(222); hold on
hsh(1)=shadedErrorBar(timeBslDown,grandMeansCorrDown(3,1:numel(timeBslDown))  , grandSemsCorrDown(1,1:numel(timeBslDown)) ,'b' ,1);
hsh(2)=shadedErrorBar(timePostDown,grandMeansCorrDown(3,numel(timeBslDown)+1:end),grandSemsCorrDown(1,numel(timeBslDown)+1:end) ,'b' ,1);
hsh(3)=shadedErrorBar(timeBslDown,grandMeansCorrDown(4,1:numel(timeBslDown))  , grandSemsCorrDown(2,1:numel(timeBslDown)) ,'r',1);
hsh(4)=shadedErrorBar(timePostDown,grandMeansCorrDown(4,numel(timeBslDown)+1:end),grandSemsCorrDown(2,numel(timeBslDown)+1:end),'r',1);


hs(3)=subplot(223); hold on
hsh(1)=shadedErrorBar(timeBslDown,grandMeansCorrDown(5,1:numel(timeBslDown))  , grandSemsCorrDown(5,1:numel(timeBslDown)) ,'b' ,1);
hsh(2)=shadedErrorBar(timePostDown,grandMeansCorrDown(5,numel(timeBslDown)+1:end),grandSemsCorrDown(5,numel(timeBslDown)+1:end) ,'b' ,1);
hsh(3)=shadedErrorBar(timeBslDown,grandMeansCorrDown(6,1:numel(timeBslDown))  , grandSemsCorrDown(6,1:numel(timeBslDown)) ,'r',1);
hsh(4)=shadedErrorBar(timePostDown,grandMeansCorrDown(6,numel(timeBslDown)+1:end),grandSemsCorrDown(6,numel(timeBslDown)+1:end),'r',1);
ylabel('10 Hz Power');

hs(4)=subplot(224); hold on
hsh(1)=shadedErrorBar(timeBslDown,grandMeansCorrDown(7,1:numel(timeBslDown))  , grandSemsCorrDown(5,1:numel(timeBslDown)) ,'b' ,1);
hsh(2)=shadedErrorBar(timePostDown,grandMeansCorrDown(7,numel(timeBslDown)+1:end),grandSemsCorrDown(5,numel(timeBslDown)+1:end) ,'b' ,1);
hsh(3)=shadedErrorBar(timeBslDown,grandMeansCorrDown(8,1:numel(timeBslDown))  , grandSemsCorrDown(6,1:numel(timeBslDown)) ,'r',1);
hsh(4)=shadedErrorBar(timePostDown,grandMeansCorrDown(8,numel(timeBslDown)+1:end),grandSemsCorrDown(6,numel(timeBslDown)+1:end),'r',1);


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
print -dpng ../figures/mean_envelopes_allRats_


%%
% % show subject averaged downsampled curve
% conds2show=[1 2 5 6];
% figure; 
% for c=1:numel(conds2show)
% hs(c)=subplot(2,2,c); hold on
% plot(timeBslDown,grandMeansCorrDown(conds2show(c),1:numel(timeBslDown)))
% plot(timePostDown,grandMeansCorrDown(conds2show(c),numel(timeBslDown)+1:end))
% if c<3, htit(c)=title(titles{conds2show(c)});  end
% %xlim(xl);
% end
% set(get(hs(1),'Ylabel'),'String','40 Hz Power');
% set(get(hs(2),'Ylabel'),'String','10 Hz Power');
% %set(hs(1:4),'Ylim',[0 70]);
% %set(hs(5:8),'Ylim',[0 200]);
% %set(hs(1:8),'Xlim',[-10 30]);
% %set(hs(1:8),'Xlim',xl);
% print -dpng ../figures/mean_envelopes_allRats_


% %%
% figure; 
% for c=1:8
% hs(c)=subplot(2,4,c); hold on
% plot(timeBsl,grandMeansCorr(c,1:numel(timeBsl)))
% plot(timePost,grandMeansCorr(c,numel(timeBsl)+1:end))
% if c<5, htit(c)=title(titles{c});  end
% xlim(xl);
% end
% set(get(hs(1),'Ylabel'),'String','40 Hz Power');
% set(get(hs(5),'Ylabel'),'String','10 Hz Power');
% set(hs(1:4),'Ylim',[0 70]);
% set(hs(5:8),'Ylim',[0 200]);
% %set(hs(1:8),'Xlim',[-10 30]);
% set(hs(1:8),'Xlim',xl);
% %print -dpng ../figures/envelopes_allRats_
% 
% 
% 
% %%
% % baseline corrected, active - sham, individual rats
% figure;
% for r=1:nRats
%     for c=1:4
%         hs2(c)=subplot(nRats,4,(r-1)*4+c); hold on
%         plot(timeBsl,   squeeze(allDataCorr(c*2-1,r,1:numel(timeBsl)))  - squeeze(allDataCorr(c*2,r,1:numel(timeBsl)))     )
%         plot(timePost,  squeeze(allDataCorr(c*2-1,r,numel(timeBsl)+1:end))  - squeeze(allDataCorr(c*2,r,numel(timeBsl)+1:end)))
%     end
% end
% 
% %%
% del1=squeeze(allDataCorr(1,:,numel(timeBsl)+1:end))  - squeeze(allDataCorr(2,:,numel(timeBsl)+1:end));
% del2=squeeze(allDataCorr(3,:,numel(timeBsl)+1:end))  - squeeze(allDataCorr(4,:,numel(timeBsl)+1:end));
% del3=squeeze(allDataCorr(5,:,numel(timeBsl)+1:end))  - squeeze(allDataCorr(6,:,numel(timeBsl)+1:end));
% del4=squeeze(allDataCorr(7,:,numel(timeBsl)+1:end))  - squeeze(allDataCorr(8,:,numel(timeBsl)+1:end));
% 
% nMinutes=size(del1,2)/fsd/60;
% for m=1:nMinutes
%     tinds=(m-1)*fsd*60+1:m*fsd*60;
%     ts1(:,m)=mean(del1(:,tinds),2);
%     ts2(:,m)=mean(del2(:,tinds),2);
%     ts3(:,m)=mean(del3(:,tinds),2);
%     ts4(:,m)=mean(del4(:,tinds),2);
%     [p1]=signrank(mean(del1(:,tinds),2));
%     [p2]=signrank(mean(del2(:,tinds),2));
%     [p3]=signrank(mean(del3(:,tinds),2));
%     [p4]=signrank(mean(del4(:,tinds),2));
%     pvals(:,m)=[p1;p2;p3;p4];
% end
%     
% 
% % set(get(hs2(1),'Title'),'String','Active 40 - Sham 40');
% % set(get(hs2(2),'Title'),'String','Active 10 - Sham 10');
% % set(get(hs2(1),'Ylabel'),'String','40 Hz Power');
% % set(get(hs2(3),'Ylabel'),'String','10 Hz Power');
% % set(hs2(1:2),'Ylim',[-40 60]);
% % set(hs2(3:4),'Ylim',[-150 150]);
% % set(hs2(1:4),'Xlim',[-10 40]);
% %print -dpng ../figures/corrected_activeminussham_envelopes_allRats_



% %%
% % baseline corrected, active - sham
% figure; 
% for c=1:4
% hs2(c)=subplot(2,2,c); hold on
% plot(timeBsl,   grandMeansCorr(c*2-1,1:numel(timeBsl))  - grandMeansCorr(c*2,1:numel(timeBsl))     )
% plot(timePost,  grandMeansCorr(c*2-1,numel(timeBsl)+1:end)  - grandMeansCorr(c*2,numel(timeBsl)+1:end))
% end
% set(get(hs2(1),'Title'),'String','Active 40 - Sham 40');
% set(get(hs2(2),'Title'),'String','Active 10 - Sham 10');
% set(get(hs2(1),'Ylabel'),'String','40 Hz Power');
% set(get(hs2(3),'Ylabel'),'String','10 Hz Power');
% set(hs2(1:2),'Ylim',[-40 60]);
% set(hs2(3:4),'Ylim',[-150 150]);
% set(hs2(1:4),'Xlim',[-10 40]);
% %print -dpng ../figures/corrected_activeminussham_envelopes_allRats_






