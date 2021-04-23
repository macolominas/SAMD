function ampOut = filterAmplt(ampIn,thre)
% Filtering in the frequency domain
%
% By Haizhao Yang

[num N] = size(ampIn);
ampOut = zeros(num,N);
for cnt = 1:num
    temp = fftshift(fft(ifftshift(ampIn(cnt,:))))/N;
    temp(min(round(N/2+thre),N):end) = 0;
    temp(1:max(1,round(N/2-thre))) = 0;
    ampOut(cnt,:) = real(fftshift(ifft(ifftshift(temp)))*N);
    ampOut(cnt,:) = ampOut(cnt,:)*norm(ampIn(cnt,:))/norm(ampOut(cnt,:));
end
      