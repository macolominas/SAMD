% Try fetas ECG

if (1)
   % close all;
    clear all;
    %% generate signal
    close all ; clear all ;
    
    [header, recorddata] = edfread('r01.edf') ;
    
    
    DS = 1;
    
    basicTF.win = 1024; %4096;
    basicTF.hop = 20; %441;
    basicTF.fs = 1000/DS;
    basicTF.fr = 0.02;
    basicTF.feat = 'SST11';
    advTF.num_tap = 1;%num_tap(cc);
    advTF.win_type = 'Gauss'; %{'Gaussian','Thomson','multipeak','SWCE'}; %}
    advTF.Smo = 1;
    advTF.Rej = 0;
    advTF.ths = 1E-9;
    advTF.HighFreq = 10/100;
    advTF.LowFreq = 0.1/100;
    advTF.lpc = 0;
    cepR.g = 0.3; % g(cc);%0.06;
    cepR.Tc=0;
    P.num_s = 1;
    P.num_c = 1;
    
    
    time_stamp = basicTF.hop/basicTF.fs;
    numGroup = 2;
    
    x0 = recorddata(2,:) ; x = zeros(size(x0)) ;
    for ii = 1:length(x)
        idx = max(1, ii-35):min(length(x), ii+35) ;
        x(ii) = x0(ii) - median(x0(idx)) ;
    end
    x0 = x ; x = resample(x0,1,DS);
%     [tfrQ4, ceps4, tceps4, tfrrQ4, rtfrQ4, tfrsqQ4, tfrtic4] = CFPH(x, basicTF, advTF, cepR, P);
%     TIME = 0:time_stamp:time_stamp*(size(rtfrQ4,2)-1) ;
%     maxT = max(TIME);
%     TIME = TIME/maxT;
%     FREQ = tfrtic4*basicTF.fs*maxT;
%     [fridge1,~,lridge] = tfridge(rtfrQ4(1:90,:),FREQ(1:90,1),0.1,'NumRidges',1,'NumFrequencyBins',10);
%     figure;plot(TIME,fridge1,'linewidth',1)
%     
%     [fridge2,~,lridge] = tfridge(rtfrQ4(80:150,:),FREQ(80:150,1),0.1,'NumRidges',1,'NumFrequencyBins',5);
%     
%     figure;plot(TIME,fridge2,'linewidth',1)
    
    sig = x(1:2^16);
    N = 2^floor(log2(length(sig)));
    sig = sig(1:N);
    x = (0:N-1)/N;
    sig = sig - mean(sig);
    
    numGroup = 1;
    opt.eps = 1e-3;
    opt.res = 1;
    opt.freq_range = [0 400];
    opt.NG = N;
    opt.dt = 1/N;
    opt.t_sc = 0.75;
    opt.NM = 0;
    opt.st = round([200]);
    opt.ed = round([300]);
    opt.num_select = numGroup;
    opt.red = 8;
    opt.C = 50;
    opt.rad = 1;
    opt.show = 1;
    [insFreq,~,insPhase,comp_select] = insInfo(sig,opt);
    opt.freq_range = [0 500];
    opt.NG = N;
    opt.dt = 1/N;
    opt.t_sc = 0.6;
    opt.NM = 0;
    opt.st = round([340]);
    opt.ed = round([380]);
    res = sig - comp_select;
    [insFreq2,~,insPhase2,comp_select2] = insInfo(res,opt);

    numGroup = 2;
    insFreq = [insFreq;insFreq2]; insPhase = [insPhase;insPhase2];
    comp_select = [comp_select;comp_select2];
    %sig = x;
    N = length(sig);
    % find peaks
    peaks = zeros(size(comp_select));
    for cnt = 1:numGroup
        peaks(cnt,:) = peakDetection(comp_select(cnt,:),insFreq(cnt,:));
    end
    
    % correct phases
     insPhase = phaseShift(insPhase,peaks);
    
    opt.maxiter = 200;
    opt.eps_error = 1e-6;
    opt.show = 0;
    opt.iterStyle = 'GS';
    opt.shapeMethod = 1;
    opt.eps_diff = 1e-6;
    opt.ampErrBandWidth = 20;
    opt.numSweep = 1;
    
    switch opt.shapeMethod
        case 1
            opt.para.Ls=1000;
            opt.para.bandWidth = 20;
            opt.para.diffeoMethod = 'nufft';
        case 2
            opt.para.nknots = 20;
            opt.para.knotremoval_factor= 1.0001;
            opt.para.order = 3;
            opt.para.Ls = 1000;
    end
    numGroup = 2;
    % test example: two components
    [shape,comp,Hcoef] = DeCom_MMD(sig,x,numGroup,ones(size(insFreq)),insFreq,insPhase,opt);
    
    save('./results/RDSA_fig10.mat','-v7.3');
end

if (1)
    load ./results/RDSA_fig10.mat;
    Ntotal = length(sig);
    N1 = N/16+1; N2 = 3*N/16;
    
    xvec = (N1:N2)/Ntotal*480;
    
    pic = figure;
    plot(insFreq(1,:),'LineWidth',2);
    hold on;
    plot(insFreq(2,:),'LineWidth',2);
    hold off;
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,0,1000]);
    title('remove'); legend('Respiratory','Cardiac');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = sprintf('./results/RDSA_fig23_%d',8+cnt-1); saveas(pic,str);
    print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
    
    %% Modified to regenerate the signal to test if it is actually accurate enough
    Sig = comp{1} + comp{2};
    
    st = 30; ed = 200;
    res = sig(N1:N2)-Sig(N1:N2);
    L = numel(res);
    resh = fftshift(fft(res));
    resh(L/2-st:L/2+st) = 0;
    %resh(1:L/2-ed) = 0; resh(L/2+ed:end) = 0;
    res = ifft(ifftshift(resh));
    res = abs(res - mean(res));
    res = wden(res,'rigrsure','h','one',8,'sym4');
    windowWidth = int16(5);
    halfWidth = windowWidth / 2;
    sz = 500;
    gaussFilter = gausswin(sz);
    gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.
    
    % Do the blur.
    size(res)
    res = conv(res, gaussFilter);
    size(res)
    res = res(sz/2:end-sz/2);
    posp = find(res>0.2);
    
    
    pic = figure;
    plot(xvec,comp{1}(N1:N2),'b'); hold on; sigy = comp{1}(N1:N2);
    axis tight; xlabel('time');ylabel('signal intensity');
    title('comp1');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-60,20]);
    xlabel('time'); ylabel('frequency');
    saveas(pic,'./results/RDSA_fig23_comp1.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig23_comp1';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(xvec,comp{2}(N1:N2),'b'); hold on; sigy = comp{2}(N1:N2);
    axis tight; xlabel('time');ylabel('signal intensity');
    title('comp2');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-60,20]);
    xlabel('time'); ylabel('intensity');
    saveas(pic,'./results/RDSA_fig23_comp2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig23_comp2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    

    
    pic = figure;
    plot(xvec,Sig(N1:N2),'r');
    pbaspect([10 1 1]);set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-60,20]);
    hold on
    plot(xvec,sig(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    saveas(pic,'./results/RDSA_fig23_sum.fig');
    title('compare reconstruction');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig23_sum';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(xvec,sig(N1:N2),'b'); hold on; sigy = sig(N1:N2);
    axis tight; xlabel('time');ylabel('signal intensity');
    title('original signal');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-60,20]);
    saveas(pic,'./results/RDSA_fig23_org.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig23_org';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(xvec,sig(N1:N2)-Sig(N1:N2),'b'); hold on; sigy = sig(N1:N2)-Sig(N1:N2);
    axis tight; xlabel('time');ylabel('signal intensity');
    title('res');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-60,20]);
    saveas(pic,'./results/RDSA_fig23_res.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig23_res';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    res = sig-Sig;
    XC = xcorr(res)/length(res);
    pic = figure;
    plot(XC,'b'); axis tight; xlabel('frequency');ylabel('power spectral density');
    saveas(pic,'./results/RDSA_fig23_psd.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig23_psd';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
end

opt.maxiter = 200;
opt.eps_error = 1e-6;
opt.show = 0;
opt.iterStyle = 'GS';
opt.shapeMethod = 1;
opt.eps_diff = 1e-6;
opt.ampErrBandWidth = 0;
opt.numSweep = 1;

% test example: two components
[shape,comp] = DeCom_MMD(sig,x,numGroup,ones(size(insFreq)),insFreq,insPhase,opt);

%% Modified to regenerate the signal to test if it is actually accurate enough
Sig = comp{1} + comp{2};

res = sig-Sig;
XC = xcorr(res)/length(res);
pic = figure;
plot(XC,'b'); axis tight; xlabel('frequency');ylabel('power spectral density');
saveas(pic,'./results/RDSA_fig23_psd2.fig');
set(gca, 'FontSize', 16);
b=get(gca);
set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
str = './results/RDSA_fig23_psd2';
print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);



XC = xcorr(randn(size(res)))/length(res);
pic = figure;
plot(XC,'b'); axis tight; xlabel('frequency');ylabel('power spectral density');
saveas(pic,'./results/RDSA_fig23_psd3.fig');
set(gca, 'FontSize', 16);
b=get(gca);
set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
str = './results/RDSA_fig23_psd3';
print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);


