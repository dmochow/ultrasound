% read output of kilosort spike sorting
clear all; close all; clc
addpath(genpath('D:\Kilosort2-master'));
addpath(genpath('D:\npy-matlab-master'));

% parameters
minfr=10;
maxfr=100;
nBins=10;

% constants
fs=30000;
tOn=2*60*fs;
tOff=5*60*fs;
tEnd=10*60*fs;

%%
ratStrs={'r027','r028','r029','r031','r032','r033','r034','r035','r036','r037',...
    'r038','r039','r040','r041','r042','r043','r044','r045','r046','r047'};
nRats=numel(ratStrs);
intStrs={'2.5Wcm2','5.0Wcm2','10.0Wcm2'};
nInts=numel(intStrs);
condStrs={'Active','Sham'};
nConds=numel(condStrs);
outDataRoot='D:\Kilosort2-results\data\Study3\raw\mat\';
rezFilename='rez2.mat';
spikeTimes=cell(nRats,nInts,nConds);
allLambdas=cell(nInts,nConds);
spikeDuration=cell(nInts,nConds);

for i=1:nInts
    for c=1:nConds
        for r=1:nRats
            
            ratStr=ratStrs{r};
            intStr=intStrs{i};
            condStr=condStrs{c};
            
            [i,c,r]
            pathToRezFile = [outDataRoot '\' ratStr '\' intStr '\' condStr '\' ];
            
            if exist(fullfile(pathToRezFile,rezFilename)) % some don't have any spikes
                
                load(fullfile(pathToRezFile,rezFilename));
                spikeTable=rez.st3;
                nClusts=max(spikeTable(:,2));
                
                tspikeTimes=cell(nClusts,1);
                lambda=zeros(nClusts,3);
                
                % fill the array
                for cl=1:nClusts
                    if 1%rez.good(cl)
                        tspikeTimes{cl,1}=spikeTable(spikeTable(:,2)==cl,1);
                    else
                        tspikeTimes{cl,1}=[];
                    end
                end
                
                % count pre/during/post
                lambda(:,1)=cellfun(@(x)sum(x<tOn),tspikeTimes)/(tOn/fs);
                lambda(:,2)=cellfun(@(x)sum(x>tOn&x<tOff),tspikeTimes)/((tOff-tOn)/fs);
                lambda(:,3)=cellfun(@(x)sum(x>tOff),tspikeTimes)/((tEnd-tOff)/fs);
                
                allLambdas{i,c}=cat(1,allLambdas{i,c},lambda);
                
                % look at templates here
                npyFilename='templates.npy';
                spikeTemplates = readNPY(fullfile(pathToRezFile,npyFilename));
                clear tspikeDuration 

                for cl=1:nClusts
                    
%                     figure(1)
                    tspikeTemplates=squeeze(spikeTemplates(cl,:,:));
                    nGoodChan=size(tspikeTemplates,2);
                    nRows=3; nCols=ceil(nGoodChan/nRows);
                    for ch=1:nGoodChan
                        twave(:,ch)=tspikeTemplates(:,ch);
                        if trapz(twave(:,ch))>0
                            twave(:,ch)=NaN*twave(:,ch);
                            tspikeDuration(cl,ch)=NaN;
                        else
                            twave(:,ch)=smoothdata(twave(:,ch));
                            [~,trough]=min(twave(:,ch));
                            twave(1:trough-1,ch)=-inf; % to ensure max comes after trough
                            [~,peak]=max(twave(:,ch));
                            tspikeDuration(cl,ch)=peak-trough;
                        end
%                         subplot(nRows,nCols,ch);
%                         plot(twave(:,ch));
                        

                       
                    end
                     
%                     for ch=1:size(spikeTemplates,3)
%                         twave=tspikeTemplates(:,ch);    
%                         if trapz(twave)>0
%                             tspikeDuration(cl,ch)=NaN;
%                         else
%                             %smtwave=smoothdata(twave);
%                             smtwave=twave;
%                             [~,trough]=min(smtwave);
%                             %smtwave(1:trough)=-inf; % to ensure max comes after trough
%                             [~,peak]=max(smtwave);
%                             tmp=peak-trough;
%                             if tmp>0
%                                 tspikeDuration(cl,ch)=tmp;
%                             else
%                                 tspikeDuration(cl,ch)=NaN;
%                             end
%                            
%                         end                  
%                     end
                end
                spikeDuration{i,c}= cat(1,spikeDuration{i,c},nanmedian(tspikeDuration,2));
                
%                 plot([trough trough],[yl(1) yl(2)],'--b')
%                 plot([peak peak],[yl(1) yl(2)],'--r')
                
            end
            
        end
    end
end


%%
% mus=cellfun(@(x)(nanmean(x,1)),allLambdas,'UniformOutput',false);
% musAcute=cellfun(@(x)(nanmean( ( x(:,2)-x(:,1))./x(:,1)) ),allLambdas,'UniformOutput',false);
% musChronic=cellfun(@(x)(nanmean( ( x(:,3)-x(:,1))./x(:,1)) ),allLambdas,'UniformOutput',false);


%%
for i=1:3
    for c=1:2
        badIdx=find(allLambdas{i,c}(:,1)<minfr | allLambdas{i,c}(:,1)>maxfr);
        frsPre{i,c}=allLambdas{i,c}(setdiff(1:size(allLambdas{i,c},1),badIdx),1);
        frsDur{i,c}=allLambdas{i,c}(setdiff(1:size(allLambdas{i,c},1),badIdx),2);
        frsPost{i,c}=allLambdas{i,c}(setdiff(1:size(allLambdas{i,c},1),badIdx),3);
    end
end

%figure
for i=1:3
    delAcuteActive{i,1}= (frsDur{i,1} - frsPre{i,1}) ./  frsPre{i,1};
    delAcuteSham{i,1}= (frsDur{i,2} - frsPre{i,2}) ./  frsPre{i,2};
    
    delChronicActive{i,1}= (frsPost{i,1} - frsPre{i,1}) ./  frsPre{i,1};
    delChronicSham{i,1}= (frsPost{i,2} - frsPre{i,2}) ./  frsPre{i,2};
    
    pvalAcute(i)=ranksum(delAcuteActive{i},delAcuteSham{i});
    pvalChronic(i)=ranksum(delChronicActive{i},delChronicSham{i});
end

%%
figure;
for i=1:3
    subplot(2,3,i); hold on
    [Na,Xa]=hist(delAcuteActive{i}*100,nBins);
    [Ns,Xs]=hist(delAcuteSham{i}*100,nBins);
    set(gca,'YScale','log');
    bar(Xa,Na,'b','FaceAlpha',0.5);
    bar(Xs,Ns,'r','FaceAlpha',0.5);
    xlabel('Percent change');
    ylabel('# of neurons')
    
    subplot(2,3,i+3); hold on
    [Na,Xa]=hist(delChronicActive{i}*100,nBins);
    [Ns,Xs]=hist(delChronicSham{i}*100,nBins);
    set(gca,'YScale','log');
    bar(Xa,Na,'b','FaceAlpha',0.5);
    bar(Xs,Ns,'r','FaceAlpha',0.5);
    xlabel('Percent change');
    ylabel('# of neurons')
end

