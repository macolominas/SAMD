function R = DeCom_whiteness(x,m)
    N = length(x);
    XC = xcorr(x)/N;
    figure;subplot(1,2,1);plot(XC);subplot(1,2,2);plot(abs(fftshift(fft(XC))/N).^2)
    XC = XC(N:end);
    lags = (1:m)+1;
    R = N*sum(XC(lags).^2)/XC(1).^2;
end