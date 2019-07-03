% 06.06.19
% adapt to handle data folders with (r021-) and without stim (r012-r019)

clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));

%%
rStr='r024';
freq=40;
condStr='active';
freqStr=[num2str(freq) 'HZ'];
subfolderStr=[condStr num2str(freq)]; 
dataFolder=['../data/' rStr '/' subfolderStr '/' ];
saveDataFileStr=[rStr '_' freqStr '_' condStr '_hilbert'];

%%
rerefChannel=NaN;
fs=30000;
dsr=60; %downsampling ratio (fsr = fs/dsr)
flHz=1; % lowest frequency that we're interested in
fhHz=250; % highest frequency we interested in
fsr=fs/dsr;
butterOrder=2;
probeMap=[17,16,18,15,19,14,20,13,...
    21,12,22,11,23,10,24,9,...
    25,8,26,7,27,6,28,5,...
    29,4,30,3,31,2,32,1]';
filenames=dir([dataFolder '*.mat']);
nFiles=numel(filenames);
for f=1:nFiles
    tFilename=upper(filenames(f).name);
    if contains(tFilename,'BASELINE') || contains(tFilename,'PRE') 
        bslFilename=tFilename;
    elseif contains(tFilename,'POSTTUS_1') || contains(tFilename,'POS_1') || contains(tFilename,'POST1')  
        post1Filename=tFilename;
    elseif contains(tFilename,'POSTTUS_2') || contains(tFilename,'POS_2') || contains(tFilename,'POST2')
        post2Filename=tFilename;
    elseif contains(tFilename,'POSTTUS_3') || contains(tFilename,'POS_3') || contains(tFilename,'POST3')
        post3Filename=tFilename;
    elseif contains(tFilename,'STIM') 
       stimFilename=tFilename;
       stimExists=1;
    else
    end
end
matFilenames{1}=fullfile(dataFolder,bslFilename);
matFilenames{2}=fullfile(dataFolder,post1Filename);
matFilenames{3}=fullfile(dataFolder,post2Filename);
matFilenames{4}=fullfile(dataFolder,post3Filename);
if stimExists
    matFilenames{5}=fullfile(dataFolder,stimFilename);
    samples={60*0*fs+1:60*10*fs; % 10 min
    60*0*fs+1:60*10*fs;
    60*0*fs+1:60*10*fs;
    60*0*fs+1:60*10*fs;
    60*0*fs+1:60*3*fs} % 3 min}; 
else
    samples={60*0*fs+1:60*10*fs; % 10 min
    60*0*fs+1:60*10*fs;
    60*0*fs+1:60*10*fs;
    60*0*fs+1:60*10*fs} % 3 min}; 
end

%%
for f=1:4 % pre + post1 + post2 + post3
    f
    load(matFilenames{f});
    if size(amplifier_data,2) < samples{f}(end)
        fprintf('Warning: missing data \n'); 
        nMissing= samples{f}(end) - size(amplifier_data,2);
        amplifier_data=cat(2,amplifier_data,zeros(size(amplifier_data,1),nMissing));
    end
    data(:,:,f)=amplifier_data(:,samples{f}); 
    clear amplifier_data;
end

%%
% bring in the data during stimulation
if stimExists
    load(matFilenames{5});
    stimData=amplifier_data(:,samples{5}); 
    clear amplifier_data;
end

%%
% rereference
if ~isnan(rerefChannel)
    data = data - repmat (data(rerefChannel,:,:),[size(data,1) 1 1 ]);
    if stimExists
        stimData = stimData - repmat (stimData(rerefChannel,:),[size(stimData,1) 1]);
    end
end

%%
% initial bandpass filter to bring data to 1-250 Hz
fl=flHz/(fs/2);
fh=fhHz/(fs/2);
[b,a]=butter(butterOrder,[fl fh]);
data=filter(b,a,data,[],2);
if stimExists
    stimData=filter(b,a,stimData,[],2);
end


%%
% downsample data (lowpass filtering already done so we can just pick every dsrth sample)
dataDown=data(:,1:dsr:end,:);
if stimExists
    stimDataDown=stimData(:,1:dsr:end);
end


%%
% free memory
clear data
clear stimData

%%
% notch filters
nHarm=floor(fhHz/60);
for h=1:nHarm
    h
    [bn,an]=butter(2,[h*60-2 h*60+2]/fsr*2,'stop');
    dataDown=filter(bn,an,dataDown,[],2);
    if stimExists
        stimDataDown=filter(bn,an,stimDataDown,[],2);
    end
end

%% 
% robust pca
for b=1:size(dataDown,3)
    tmp=squeeze(dataDown(:,:,b)).'; % time x space
    [A_hat E_hat iter] = inexact_alm_rpca(tmp);
    dataDown(:,:,b)=A_hat.';
end

if stimExists
    tmp=stimDataDown.'; % time x space
    [A_hat E_hat iter] = inexact_alm_rpca(tmp);
    stimDataDown=A_hat.';
end

%%
% tmpPre=dataDown(:,:,1); tmpPost=dataDown(:,:,2:nFilenames-1);
% allDataDown= [ tmpPre(:,:) stimDataDown tmpPost(:,:) ]; % combined 43 minute block
% 

%%
% bandpass, hilbert, abs, save

% design filters
[b40,a40]=butter(2,[36 44]/fsr*2);
if ~isstable(b40,a40), error('filter not stable bro'); end
[b10,a10]=butter(2,[6 14]/fsr*2);
if ~isstable(b10,a10), error('filter not stable bro'); end

%%
% apply filters
dataDown40=filter(b40,a40,dataDown,[],2);
dataDown10=filter(b10,a10,dataDown,[],2);

% tmp=(hilbert(allDataDown40.')).';
% env40=abs(tmp);
% 
% tmp=(hilbert(allDataDown10.')).';
% env10=abs(tmp);

for bl=1:size(dataDown,3)
    tmp=(hilbert(dataDown40(:,:,bl).')).';
    env40(:,:,bl)=abs(tmp);
    
    tmp=(hilbert(dataDown10(:,:,bl).')).';
    env10(:,:,bl)=abs(tmp);
end

if stimExists
    stimDataDown40=filter(b40,a40,stimDataDown,[],2);
    tmp=(hilbert(stimDataDown40.')).';
    stimenv40=abs(tmp);
    
    stimDataDown10=filter(b10,a10,stimDataDown,[],2);
    tmp=(hilbert(stimDataDown10.')).';
    stimenv10=abs(tmp);
end


%%
% subsample envelope to 100 Hz sampling
% THIS ASSUMES FS=30000 and FSR=500
% FIXTHIS
dsr2=10; 
fsr2=fsr/dsr2;
[bd,ad]=butter(2,50/(fsr*2),'low');
if ~isstable(bd,ad), error('filter not stable bro'); end
env40d=filter(bd,ad,env40,[],2);
env10d=filter(bd,ad,env10,[],2);
env40ds=env40d(:,1:dsr2:end,:);
env10ds=env10d(:,1:dsr2:end,:);

if stimExists
    stimenv40d=filter(bd,ad,stimenv40,[],2);
    stimenv10d=filter(bd,ad,stimenv10,[],2);
    stimenv40ds=stimenv40d(:,1:dsr2:end);
    stimenv10ds=stimenv10d(:,1:dsr2:end);
end


%%
% save for input to a second script
if stimExists
    save(fullfile(dataFolder,[date '-' saveDataFileStr]),'fsr','fsr2',...
    'env40ds','env10ds','stimenv40ds','stimenv10ds','rerefChannel','butterOrder','flHz','fhHz','probeMap');
else
    save(fullfile(dataFolder,[date '-' saveDataFileStr]),'fsr','fsr2',...
    'env40ds','env10ds','rerefChannel','butterOrder','flHz','fhHz','probeMap');
end
    

