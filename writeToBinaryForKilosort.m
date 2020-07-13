% convert smartbox files to binary for kilosort
%clear all; close all; clc
%addpath(genpath('D:\Kilosort2-master'));
%addpath(genpath('D:\npy-matlab-master'));

function writeToBinaryForKilosort(pathToReadData,pathToWriteData,ratStr,condStr)

%%
% bring in channel map so that data is ordered from bottom to top
probeMap=[17,16,18,15,19,14,20,13,...
    21,12,22,11,23,10,24,9,...
    25,8,26,7,27,6,28,5,...
    29,4,30,3,31,2,32,1]';
probeMapReverse=probeMap(end:-1:1);
%pathToData='Z:\ULTRASOUND\data\Study3\raw\mat\r027\2.5Wcm2';
filenames=dir(fullfile(pathToReadData,[ratStr '_' condStr '*.mat']));
filename=filenames(1).name; %  it better be the first one!
%filename='r027_Active_20190826_174100.mat';
vars=load(fullfile(pathToReadData,filename));
fs=30000;
%%
data=vars.amplifier_data(probeMapReverse,1:10*60*fs);
dataI=int16(data);
binFilename=cat(2,filename(1:end-4),'_reordered.bin');
fid = fopen(fullfile(pathToWriteData,binFilename), 'w'); 
fwrite(fid, dataI, 'int16');
fclose(fid);
