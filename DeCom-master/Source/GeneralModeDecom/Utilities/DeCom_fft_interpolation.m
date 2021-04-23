function y = DeCom_fft_interpolation(sig,N)
% Interpolation via padding in the frequency domain using FFT
%
% By Haizhao Yang

[m,n] = size(sig);
N0 = max(m,n);
if N0<N
    sigh = fftshift(fft(sig));
    if m > n
        sighNew = zeros(N,1);
    else
        sighNew = zeros(1,N);
    end
    st = N/2+1-floor(N0/2);
    sighNew(st:st+N0-1) = sigh;
    y = real(ifft(ifftshift(sighNew))*(sqrt(N))/sqrt(N0))*2;
else
    y = sig;
end