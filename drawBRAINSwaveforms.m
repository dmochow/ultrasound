clear all; close all; clc

% draw conventional and AM US

%% conventional
fsx=0.3/20; nPulses=10;
sp=cos(2*pi*fsx*(1:100));
spp=[sp zeros(size(sp))];
sppr=repmat(spp,[1 nPulses]);
figure
hs(1)=subplot(231);
plot(sppr); axis off
title('Stimulus 1','FontWeight','normal');
text(2150,0,'ISI','FontSize',14);
hs(2)=subplot(232);
plot(sppr); axis off
text(2250,0.1,'...','FontSize',20);
title('Stimulus 2','FontWeight','normal');
hs(3)=subplot(233);
plot(sppr); axis off
title('Stimulus N','FontWeight','normal');

pos=get(hs(1),'Position');
set(hs(1),'Position',[0.025 pos(2) 0.25 0.2]);

pos=get(hs(2),'Position');
set(hs(2),'Position',[0.35 pos(2) 0.25 0.2]);

pos=get(hs(3),'Position');
set(hs(3),'Position',[0.7 pos(2) 0.25 0.2]);

%% proposed
fsx=1/20; 
fsm=0.05/20;
s=abs(sin(2*pi*fsm*(1:1400))).*cos(2*pi*fsx*(1:1400));
sh=hilbert(s);
hs(4)=subplot(212); hold on
plot(s);
plot(abs(sh));
axis off;
title('Amplitude Modulated FUS','FontWeight','normal','FontSize',16);
hlg=legend('Fast carrier (2 MHz)','Slow envelope (40 Hz)');
set(hlg,'box','off');
set(hlg,'orientation','horizontal');
set(hlg,'location','southeast');
lgpos=get(hlg,'Position');
set(hlg,'Position',[lgpos(1)-0.15 lgpos(2)+0.15 lgpos(3) lgpos(4)]);

pos=get(hs(4),'Position');
set(hs(4),'Position',[0.025 0.3 0.95 0.2]);

htext=text(500,4.4,'Conventional FUS','FontSize',16);

%%

hs=sublabel([hs(1) hs(4)],-20,-10,'FontWeight','normal','FontSize',16);
%%
print -dpng ../figures/BRAINSwaveforms
crop('../figures/BRAINSwaveforms.png',0);



% 
% fo=0.5e6;
% cpp=100;
% nPulses=80;
% prf=1.5e3;
% fs=16*fo;
% onSamples=round(cpp/fo*fs);
% offSamples=round( (1/prf-cpp/fo)*fs );
% M=repmat( [ones(onSamples,1); zeros(offSamples,1)],80,1);
% nSamples=(onSamples+offSamples)*nPulses;
% timeMs=(0:nSamples-1)/fs*1000;
% s=cos(2*pi*fo/(4*fo)*(0:nSamples-1));
% s=s(:).*M;
% hs(1)=subplot(231);
% plot(timeMs,s); axis off
% title('FUS Stimulus 1');
% 
% hs(2)=subplot(232);
% plot(timeMs,s); axis off
% title('FUS Stimulus 2');
% 
% hs(3)=subplot(233);
% plot(timeMs,s); axis off
% title('FUS Stimulus N');
% 
% %hs(2)=subplot(222);
% %plot(timeMs,s);
% %axis tight
% %xlim([0 0.025]);
% 
% %% AM-FUS
% fo=2e6;
% fs=16*fo; 
% fm=40; % modulation frequency
% nSamples=round(20/fm*fs); % number of samples to show on the plot
% timeMs=(0:nSamples-1)/fs*1000;
% s2=abs(cos(2*pi*fm/fs*(0:nSamples-1))).*cos(2*pi*fo/fs*(0:nSamples-1));
% hs(4)=subplot(212);
% plot(timeMs,s2);
% axis off
% title('AM-FUS Stimulus');
% 
% % hs(6)=subplot(236);
% % plot(timeMs,s2);
% % axis off
% 
% %hs(4)=subplot(224);
% %plot(timeMs,s2); xlim([0 0.01]);
% 

