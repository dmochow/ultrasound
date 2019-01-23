clear all; close all; clc

load ../data/ONDA/firstHydrophoneTrace.mat

whos
VpPa=2.573e-007;
dataPa=data/1e6/VpPa;  % data in uV (?)
time=(0:numel(data)-1)/fsHardware;
[b,a]=fir1(100,(1e6/(fsHardware/2)));
filtData=filter(b,a,dataPa);
figure; %subplot(211);
plot(time*1e6,filtData/1e3);
xlabel('Time (\mus)','FontSize',16);
ylabel('Pressure (kPa)','FontSize',16);
box off
set(gca,'Xtick',[0:1:5]);
set(gca,'Ytick',[-10:5:15]);
fpos=get(gca,'Position');
set(gca,'Position',[fpos(1) fpos(2) fpos(3) fpos(3)*3/4]);
print -dpng ../figures/firstHydrophoneTrace
crop('../figures/firstHydrophoneTrace.png',0);