clear all; close all; clc

% draw conventional and AM US

%% conventional
fo=0.5e6;
cpp=100;
nPulses=80;
prf=1.5e3;
fs=16*fo;
onSamples=round(cpp/fo*fs);
offSamples=round( (1/prf-cpp/fo)*fs );
M=repmat( [ones(onSamples,1); zeros(offSamples,1)],80,1);
nSamples=(onSamples+offSamples)*nPulses;
timeMs=(0:nSamples-1)/fs*1000;
s=cos(2*pi*fo/(4*fo)*(0:nSamples-1));
s=s(:).*M;
hs(1)=subplot(231);
plot(timeMs,s); axis off
title('FUS Stimulus 1');

hs(2)=subplot(232);
plot(timeMs,s); axis off
title('FUS Stimulus 2');

hs(3)=subplot(233);
plot(timeMs,s); axis off
title('FUS Stimulus N');

%hs(2)=subplot(222);
%plot(timeMs,s);
%axis tight
%xlim([0 0.025]);

%% AM-FUS
fo=2e6;
fs=16*fo; 
fm=40; % modulation frequency
nSamples=round(20/fm*fs); % number of samples to show on the plot
timeMs=(0:nSamples-1)/fs*1000;
s2=abs(cos(2*pi*fm/fs*(0:nSamples-1))).*cos(2*pi*fo/fs*(0:nSamples-1));
hs(4)=subplot(212);
plot(timeMs,s2);
axis off
title('AM-FUS Stimulus');

% hs(6)=subplot(236);
% plot(timeMs,s2);
% axis off

%hs(4)=subplot(224);
%plot(timeMs,s2); xlim([0 0.01]);


