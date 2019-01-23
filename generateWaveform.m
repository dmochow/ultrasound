clear all; close all; clc

fs=5e6; % sampling rate
fc=0.5e6; % carrier frequency
prf=1.5e3;  % pulse-repetition-frequency
cpp=100; % number of cycles per pulse
nPulses=80; 

nSamplesPerPulse=round(fs/prf);
nSamplesPerCycle=round(fs/fc);
pulse=zeros(nSamplesPerPulse,1);
pulse(1:cpp*nSamplesPerCycle)=sin(2*pi*fc/fs*[0:cpp*nSamplesPerCycle-1]);
timep=(0:length(pulse)-1)/fs*1000;
stimulus=repmat(pulse,nPulses,1);
time=(0:length(stimulus)-1)/fs*1000;

%
nfft=2^nextpow2(nSamplesPerPulse*nPulses);
freqs=(0:nfft-1)/nfft*fs;
S=fft(stimulus,nfft);



figure; 
subplot(311); plot(timep,pulse); xlabel('Time (ms)'); title('Pulse');
subplot(312); plot(time,stimulus); xlabel('Time (ms)'); title('Stimulus');
subplot(313); plot(freqs,abs(S)); xlim([0 fs/2]); xlabel('Frequency (Hz)');


