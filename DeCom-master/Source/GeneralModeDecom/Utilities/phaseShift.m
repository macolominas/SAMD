function pOut = phaseShift(pIn,peaks)
% use weighted average to ajust input phases
%
% By Haizhao Yang

pOut = zeros(size(pIn));
[m,n] = size(pIn);
for cnt1 = 1:m % m components
    p0 = pIn(cnt1,:);
    pos = find(peaks(cnt1,:));
    pos = pos - pos(1)+1;
    numPeak = numel(pos);
    for cnt2 = 1:numPeak % need to adjust numPeak intervals
        st = pos(cnt2);
        if cnt2<numPeak
            ed = pos(cnt2+1);
        else
            ed = n;
        end
        p1 = p0(st:ed)-p0(st)+cnt2-1;
        p2 = p0(st:ed)-p0(ed)+cnt2;
        pOut(cnt1,st:ed) = avgSinWeight(p1,p2);
    end
end
end

function avg = avgSinWeight(a,b)
N = length(a);
w = sin(2*pi*(0:1/(N-1):1)).^2;
avg = w.*a+(1-w).*b;
end