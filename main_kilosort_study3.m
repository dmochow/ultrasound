% run kilosort2 on all study 3 rats
clear all; close all; clc

%%
addpath(genpath('D:\Kilosort2-master')) % path to kilosort folder
addpath('D:\npy-matlab-master') % for converting to Phy
pathToYourConfigFile = 'D:\Kilosort2-results'; % take from Github folder and put it somewhere else (together with the master_file)
chanMapFile = fullfile('Z:\ULTRASOUND\data','study3chanMap.mat');
rootH = 'D:\Kilosort2-results\tmpBinary'; % path to temporary binary file (same size as data, should be on fast SSD)


%%
myDataRoot='Z:\ULTRASOUND\data\Study3\raw\mat';
outDataRoot='D:\Kilosort2-results\data\Study3\raw\mat';
mkdir(outDataRoot);
ratStrs={'r027','r028','r029','r031','r032','r033','r034','r035','r036','r037',...
    'r038','r039','r040','r041','r042','r043','r044','r045','r046','r047'};
nRats=numel(ratStrs);
intStrs={'2.5Wcm2','5.0Wcm2','10.0Wcm2'};
nInts=numel(intStrs);
condStrs={'Active','Sham'};
nConds=numel(condStrs);

%%
% rootZ now changes dynamically
%rootZ = 'Z:\ULTRASOUND\data\Study3\raw\mat\r027\2.5Wcm2'; % the raw data binary file is in this folder

for r=4:nRats
    for i=1:nInts
        for c=1:nConds
            
            ratStr=ratStrs{r};
            intStr=intStrs{i};
            condStr=condStrs{c};
            
            inrootZ = [myDataRoot '\' ratStr '\' intStr '\' ];
            rootZ = [outDataRoot '\' ratStr '\' intStr '\' condStr '\' ];
            mkdir(rootZ);
            writeToBinaryForKilosort(inrootZ,rootZ,ratStr,condStr);
            
            ops.trange = [0 Inf]; % time range to sort
            ops.NchanTOT    = 32; % total number of channels in your recording
            
            run(fullfile(pathToYourConfigFile, 'study3Config.m'))
            ops.fproc       = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
            ops.chanMap = fullfile('Z:\ULTRASOUND\data', 'study3chanMap.mat');
            
            %% this block runs all the steps of the algorithm
            fprintf('Looking for data inside %s \n', rootZ)
            
            % find the binary file
            fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
            ops.fbinary = fullfile(rootZ, fs(1).name);
            
            % preprocess data to create temp_wh.dat
            rez = preprocessDataSub(ops);
            
            try
                % time-reordering as a function of drift
                rez = clusterSingleBatches(rez);
                
                % saving here is a good idea, because the rest can be resumed after loading rez
                save(fullfile(rootZ, 'rez.mat'), 'rez', '-v7.3');
                
                % main tracking and template matching algorithm
                rez = learnAndSolve8b(rez);
                
                % final merges
                rez = find_merges(rez, 1);
                
                % final splits by SVD
                rez = splitAllClusters(rez, 1);
                
                % final splits by amplitudes
                rez = splitAllClusters(rez, 0);
                
                % decide on cutoff
                rez = set_cutoff(rez);
                
                fprintf('found %d good units \n', sum(rez.good>0))
                
                % write to Phy
                fprintf('Saving results to Phy  \n')
                rezToPhy(rez, rootZ);
                
                %% if you want to save the results to a Matlab file...
                
                % discard features in final rez file (too slow to save)
                rez.cProj = [];
                rez.cProjPC = [];
                
                % final time sorting of spikes, for apps that use st3 directly
                [~, isort]   = sortrows(rez.st3);
                rez.st3      = rez.st3(isort, :);
                
                % save final results as rez2
                fprintf('Saving final results in rez2  \n')
                fname = fullfile(rootZ, 'rez2.mat');
                save(fname, 'rez', '-v7.3');
                
            catch
                % do nothing
            end
        end
    end
end