function peaks = peakDetection(f,freq)
% assume that the real signal f is defined in [0,1]
%
% This code detects peaks of the given signal f
%
% By Haizhao Yang

N = length(f);
peaks = zeros(1,N);

th = N/2;

rng = zeros(1,N);
for j = 1:N
    rng(j) = floor(th/freq(j));
    index = max(j-rng(j),1):min(j+rng(j),N);
    
    if(max(f(index))==f(j))
        peaks(j) = 1;
    end
end

% remove fake peaks
I = find(peaks);
d = diff(I);
num = numel(d);
L = floor(N/num);
pos = 1:L:N;
if length(pos)>length(d)
    pos = pos(1:length(d));
else if length(pos)<length(d)
        pos = [pos repmat(pos(end),[1 length(pos)-length(d)])];
    end
end
peaks(I(d<rng(pos)))=0;
peaks = fftshift(peaks);
I = find(peaks);
d = diff(I);
num = numel(d);
L = floor(N/num);
pos = 1:L:N;
if length(pos)>length(d)
    pos = pos(1:length(d));
else if length(pos)<length(d)
        pos = [pos repmat(pos(end),[1 length(pos)-length(d)])];
    end
end
peaks(I(d<rng(pos)))=0;
peaks = ifftshift(peaks);