function sigma = Chi2Estimation(y)
% estimate the standard deviation of noise from the data
%
% By Haizhao Yang

n = length(y);
yfft=fft(y(:)) / sqrt(n);
nremoved = max([1 ceil((n-1)/4)]);
nremoved = min(nremoved,floor((n-1)/2));
yfft([1+(0:nremoved-1) end-(0:nremoved-1)]) = [];
m = median([real(yfft); imag(yfft)].^2);
sigma = sqrt(m/(1-2/3+4/27-8/729));

end