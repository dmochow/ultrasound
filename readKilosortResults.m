% read output of kilosort spike sorting
clear all; close all; clc
addpath(genpath('D:\Kilosort2-master'));
addpath(genpath('D:\npy-matlab-master'));

% constants
fs=30000;
tOn=2*60*fs;
tOff=5*60*fs;
tEnd=10*60*fs;

pathToRezFile='D:\Kilosort2-results\data\Study3\raw\mat\r027\2.5Wcm2\Active';
rezFilename='rez2.mat';
load(fullfile(pathToRezFile,rezFilename));

spikeTable=rez.st3;

% initialize array holding results
nClusts=max(spikeTable(:,2));
spikeTimes=cell(nClusts,1);

% fill the array
for cl=1:nClusts
    spikeTimes{cl,1}=spikeTable(spikeTable(:,2)==cl,1);
end

% count pre/during/post
lambda=zeros(nClusts,3);
lambda(:,1)=cellfun(@(x)sum(x<tOn),spikeTimes)/(tOn/fs);
lambda(:,2)=cellfun(@(x)sum(x>tOn&x<tOff),spikeTimes)/((tOff-tOn)/fs);
lambda(:,3)=cellfun(@(x)sum(x>tOff),spikeTimes)/((tEnd-tOff)/fs);


% %%
% % 07.02.20
% % npy structures
% pathToData='Z:\ULTRASOUND\data\Study3\raw\mat\r027\2.5Wcm2';
% npyFilename='templates.npy';
% spikeTemplates = readNPY(fullfile(pathToData,npyFilename));
% 
% 
% 
%%
n=11;
Wraw = rez.Wrot'*( squeeze(rez.U(:,n,:)) * squeeze(rez.W(:,n,:))');