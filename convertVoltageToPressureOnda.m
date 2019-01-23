clear all; close all; clc

v=linspace(0.001,0.1,100);
dbre=-254.87;
tmp=10^(dbre/20); % 1uPa corresponds to this many volts
pr=(v/tmp)/(1e6*1e3);  % divide by 1e6 to get into Pa, and by 1e3 to get into kPa
figure; plot(v,pr); xlabel('Voltage (volts)'); ylabel('Pressure (kPa)'); 

tmp