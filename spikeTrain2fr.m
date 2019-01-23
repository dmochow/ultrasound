function fr = spikeTrain2fr(spikeTrain,nKernel,fs)
% convert a time series of 0s and 1s where 1s indicate a spike to a
% smoothed firing rate time series
% NB: assumes time is dimension 2

if nargin<3, fs=6000; end
if nargin<2, nKernel=30*70; end
if nargin<1, error('At least one argument needed in spikeTrain2fr'); end

sigma=0.01; % TODO: make this a parameter
%nKernel=30*70;
tKernel=(-nKernel/2+1:1:nKernel/2)/fs;
kernel=1/sqrt(2*pi*sigma.^2)*exp(-(tKernel.^2/(2*sigma^2)));
nPad=sum(tKernel<0);
% figure; plot(tKernel,kernel,'k');

%%
fr=zeros(size(spikeTrain));
for ch=1:size(spikeTrain,1)
    for bl=1:size(spikeTrain,3)
        paddedData=[zeros(1,nPad) spikeTrain(ch,:,bl)];
        smoothedData=conv(paddedData,kernel,'same');
        fr(ch,:,bl)=smoothedData(nPad+1:end);
    end
end
end

