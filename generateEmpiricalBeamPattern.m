clear all; close all; clc

pathToData='../data/watertank';
filenames={'2MHz-depth1mm-05-Nov-2018.mat'};
fo=2e6;
VpPa=1.806E-007;

for i=1:4
    if i==1
        str='';
    else
        str=['lat' num2str(i-1) 'mm'];
    end
    
    for j=1:20
        filename=['2MHz-depth' num2str(j) 'mm-' str '05-Nov-2018.mat'];
        load(fullfile(pathToData,filename));
        fs=1/mean(diff(time));
        [M,phi] = fitSineWave(scaledData,fo,fs);
        ampPa=(M/VpPa)/1000; % units are now kPa 
        Z(i,j)=ampPa;
    end
end

%%
X=0:3; Y=1:20;
figure;
pcolor(X,Y,Z'); shading interp 
xlabel('radius (mm)')
ylabel('depth (mm)')
hcb=colorbar;
set(get(hcb,'Ylabel'),'string','kPa');
%imagesc([X(1) X(end)],[Y(1) Y(end)],Z'); 
print -dpng ../figures/empiricalBeamPattern
