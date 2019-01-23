function [spikes,waveform,mua,thresh] = filterDetect(data,fs)
%UNTITLED2 denoise data and then call detect spikes
if nargin<2, fs=30000; end
if size(data,1)>size(data,2), data=data.'; warning('Transposing data'); end

% spike detection parameters
threshStds=5; % spike detection threshold (multiple of median(abs)/0.6745)
spikeDirStr='neg'; 

% specify some other 'rammies
flHz=300; % lowest frequency that we're interested in
fhHz=4000; % highest frequency we interested in
butterOrder=4;
fl=flHz/(fs/2); % trust me
fh=fhHz/(fs/2);
[b,a]=butter(butterOrder,[fl fh]);
data=filter(b,a,data,[],2);

% notch filter to remove 60 Hz + harmonics
nHarm=floor(fhHz/60); 
for h=1:nHarm
    [bn,an]=butter(2,[h*60-2 h*60+2]/fs*2,'stop');
    data=filter(bn,an,data,[],2);
end

[spikes,waveform,mua,thresh] = detectSpikes(data,threshStds,fs,spikeDirStr);
end

