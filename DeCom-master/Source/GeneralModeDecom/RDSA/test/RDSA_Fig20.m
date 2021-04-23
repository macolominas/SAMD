% Multiresolution mode decomposition for a real PPG signal.
% Figure 20-22 in the paper
% "A Fast Algorithm for Multiresolution Mode Decomposition"
% by Haizhao Yang.

if (1)
    close all;
    clear all;
    %% generate signal
    load 0009_8min.mat;
    N = 2^16;%17;
    x = (0:N-1)/N;
    sig = signal.pleth.y(1:N)';
    sig = sig - mean(sig);
    
    numGroup = 2;
    opt.eps = 1e-3;
    opt.res = 0.25;
    opt.freq_range = [0 N/128];
    opt.NG = N;
    opt.dt = 1/N;
    opt.t_sc = 0.5;
    opt.NM = 0;
    opt.st = round([50 1200 ]/16/opt.res)*N/2^14;
    opt.ed = round([400 2000 ]/16/opt.res)*N/2^14;
    opt.num_select = numGroup;
    opt.red = 8;
    opt.C = 100;
    opt.rad = 1.5;
    opt.show = 0;
    [insFreq,insAmp,insPhase,comp_select] = insInfo(sig,opt);
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
    opt.ampErrBandWidth = 40;
    opt.numSweep = 10;
    
    switch opt.shapeMethod
        case 1
            opt.para.Ls=1000;
            opt.para.bandWidth = 10;
            opt.para.diffeoMethod = 'nufft';
        case 2
            opt.para.nknots = 20;
            opt.para.knotremoval_factor= 1.0001;
            opt.para.order = 3;
            opt.para.Ls = 1000;
    end
    
    % test example: two components
    [shape,comp,Hcoef] = DeCom_MMD(sig,x,numGroup,ones(size(insAmp)),insFreq,insPhase,opt);
    Hcoef{1}.a0
    Hcoef{2}.a0
    Hcoef{1}.a
    Hcoef{2}.a
    Hcoef{1}.b
    Hcoef{2}.b
    save('./results/RDSA_fig20.mat','-v7.3');
end

if (1)
    load ./results/RDSA_fig20.mat;
    Ntotal = length(signal.pleth.y);
    N1 = N/4+1; N2 = 3*N/4;
    
    xvec = (N1:N2)/Ntotal*480;
    
    pic = figure;
    plot(xvec,insFreq(1,N1:N2),'LineWidth',2);
    hold on;
    plot(xvec,insFreq(2,N1:N2),'LineWidth',2);
    hold off;
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,0,1000]);
    title('remove'); legend('Respiratory','Cardiac');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = sprintf('./results/RDSA_fig20_%d',8+cnt-1); saveas(pic,str);
    print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
    
    
    pic = figure;
    plot(xvec,insAmp(1,N1:N2),'LineWidth',2);
    hold on;
    plot(xvec,insAmp(2,N1:N2),'LineWidth',2);
    hold off;
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,0,10]);
    title('remove'); legend('Respiratory','Cardiac');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = sprintf('./results/RDSA_fig20_%d',10+cnt-1); saveas(pic,str);
    print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
    
    
    for cnt1 = 1:numGroup
        LL = length(shape{cnt1}.s0);
        pic = figure;plot(0:1/LL:(1-1/LL),shape{cnt1}.scn{1}*Hcoef{cnt1}.a(1));
        title('remove');axis square;
        set(gca, 'FontSize', 32);
        b=get(gca);
        set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
        tit = sprintf('./results/RDSA_fig20_%d_shapeCn_%d.fig',cnt1,1);
        saveas(pic,tit);
        str = sprintf('./results/RDSA_fig20_%d_shapeCn_%d',cnt1,1);
        print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
        
        LL = length(shape{cnt1}.s0);
        pic = figure;plot(0:1/LL:(1-1/LL),shape{cnt1}.scn{2}*Hcoef{cnt1}.a(2));
        title('remove');axis square;
        set(gca, 'FontSize', 32);
        b=get(gca);
        set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
        tit = sprintf('./results/RDSA_fig20_%d_shapeCn_%d.fig',cnt1,-1);
        saveas(pic,tit);
        str = sprintf('./results/RDSA_fig20_%d_shapeCn_%d',cnt1,-1);
        print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
        
        LL = length(shape{cnt1}.s0);
        pic = figure;plot(0:1/LL:(1-1/LL),shape{cnt1}.ssn{1}*Hcoef{cnt1}.b(1));
        title('remove');axis square;
        set(gca, 'FontSize', 32);
        b=get(gca);
        set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
        tit = sprintf('./results/RDSA_fig20_%d_shapeSn_%d.fig',cnt1,1);
        saveas(pic,tit);
        str = sprintf('./results/RDSA_fig20_%d_shapeSn_%d',cnt1,1);
        print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
        
        
        LL = length(shape{cnt1}.s0);
        pic = figure;plot(0:1/LL:(1-1/LL),shape{cnt1}.ssn{2}*Hcoef{cnt1}.b(2));
        title('remove');axis square;
        set(gca, 'FontSize', 32);
        b=get(gca);
        set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
        tit = sprintf('./results/RDSA_fig20_%d_shapeSn_%d.fig',cnt1,-1);
        saveas(pic,tit);
        str = sprintf('./results/RDSA_fig20_%d_shapeSn_%d',cnt1,-1);
        print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
        
        LL = length(shape{cnt1}.s0);
        pic = figure;plot(0:1/LL:(1-1/LL),shape{cnt1}.s0*Hcoef{cnt1}.a0);
        title('remove');axis square;
        set(gca, 'FontSize', 32);
        b=get(gca);
        set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
        tit = sprintf('./results/RDSA_fig20_%d_shape_%d.fig',cnt1,0);
        saveas(pic,tit);
        str = sprintf('./results/RDSA_fig20_%d_shape_%d',cnt1,0);
        print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
    end
    
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
    gaussFilter = gausswin(sz)
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
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/RDSA_fig20_comp1.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig20_comp1';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    
    pic = figure;
    posp2 = find(posp<10000);
    mpos = min(posp(posp2)); Mpos = max(posp(posp2)); sigy = comp{1}(N1:N2);
    plot(xvec(mpos:Mpos),sigy(mpos:Mpos),'b'); hold on;
    plot(xvec(posp(posp2)),sigy(posp(posp2)),'k.');axis tight; xlabel('time');ylabel('signal intensity');
    % pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/RDSA_fig20_comp1k.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig20_comp1k';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(xvec,comp{2}(N1:N2),'b'); hold on; sigy = comp{2}(N1:N2);
    axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/RDSA_fig20_comp2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig20_comp2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    posp2 = find(10000<posp & posp<20000);
    mpos = min(posp(posp2)); Mpos = max(posp(posp2)); sigy = comp{2}(N1:N2);
    plot(xvec(mpos:Mpos),sigy(mpos:Mpos),'b'); hold on;
    plot(xvec(posp(posp2)),sigy(posp(posp2)),'k.');axis tight; xlabel('time');ylabel('signal intensity');
    % pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/RDSA_fig20_comp2k.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig20_comp2k';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(xvec,Sig(N1:N2),'r');
    pbaspect([10 1 1]);set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    hold on
    plot(xvec,sig(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    saveas(pic,'./results/RDSA_fig20_sum.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig20_sum';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(xvec,sig(N1:N2),'b'); hold on; sigy = sig(N1:N2);
    axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/RDSA_fig20_org.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig20_org';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(xvec,sig(N1:N2)-Sig(N1:N2),'b'); hold on; sigy = sig(N1:N2)-Sig(N1:N2);
    axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/RDSA_fig20_res.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig20_res';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    res = sig-Sig;
    XC = xcorr(res)/length(res);
    pic = figure;
    plot(XC,'b'); axis tight; xlabel('frequency');ylabel('power spectral density');
    saveas(pic,'./results/RDSA_fig20_psd.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig20_psd';
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
[shape,comp] = DeCom_MMD(sig,x,numGroup,insAmp,insFreq,insPhase,opt);

%% Modified to regenerate the signal to test if it is actually accurate enough
Sig = comp{1} + comp{2};

res = sig-Sig;
XC = xcorr(res)/length(res);
pic = figure;
plot(XC,'b'); axis tight; xlabel('frequency');ylabel('power spectral density');
saveas(pic,'./results/RDSA_fig20_psd2.fig');
set(gca, 'FontSize', 16);
b=get(gca);
set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
str = './results/RDSA_fig20_psd2';
print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);



XC = xcorr(randn(size(res)))/length(res);
pic = figure;
plot(XC,'b'); axis tight; xlabel('frequency');ylabel('power spectral density');
saveas(pic,'./results/RDSA_fig20_psd3.fig');
set(gca, 'FontSize', 16);
b=get(gca);
set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
str = './results/RDSA_fig20_psd3';
print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);


