% A motion artifact contaminated ECG signal to validate the effectness of
% the adaptive mode decomposition model. Fig 12-14 in the following paper.
%
% Reference:
% "Multiresolution Mode Decomposition (MMD) for Adaptive Time Series Analysis"
% by Haizhao Yang.

if (1)
    close all;
    clear all;
    %% generate signal
    %[sig,~] = plotATM('physionet_macecgdb_test01_00sm');
    [sig,~] = plotATM('physionet_macecgdb_test03_90sm');
    close all;
    N = 4000;
    sig = sig(2,1:N);
    x = (0:N-1)/N;
    
    Trend = zeros(size(sig)) ;
    for ii = 1: length(Trend)
        idx = [max(1,ii-50):min(length(Trend),ii+50)] ;
        Trend(ii) = median(sig(idx)) ;
    end
    sig = sig - Trend ;
    sig = sig - mean(sig);
    sig = sig/max(abs(sig));
    
    numGroup = 1; % number of instantaneous frequencies extracted from sigal
    opt.eps = 1e-3;
    opt.res = 0.1; % grid size in the frequency domain in the SST
    opt.freq_range = [0 N/32]; % the frequency range of instantaneous frequencies
    opt.NG = N; % number of grids in the time domain in the SST
    opt.dt = 1/N;
    opt.t_sc = 0.5;% [0.5,1]
    opt.NM = 0;
    opt.st = [8]/opt.res;
    opt.ed = [14]/opt.res;
    opt.num_select = numGroup;
    opt.red = 8; % 0.5-1.5
    opt.C = 100;
    opt.rad = 1;
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
    opt.shapeMethod = 2;
    opt.eps_diff = 1e-6;
    opt.numSweep = 10;
    
    switch opt.shapeMethod
        case 1
            opt.para.Ls=5000;
            opt.para.bandWidth = 100;
        case 2
            opt.para.nknots = 20;
            opt.para.knotremoval_factor= 1.0001;
            opt.para.order = 3;
            opt.para.Ls = 1000;
    end
    
    
    Ntotal = N;
    N1 = 1; N2 = N;
    save('./results/MMD_fig15.mat','-v7.3');
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,sig(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-5,5]);
    str = sprintf('./results/MMD_fig15_org');
    saveas(pic,[str,'.fig']);
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end

if (1)
    load ./results/MMD_fig15.mat;
    vec = [0,5,20,40];
    shape = cell(1,numel(vec));
    comp = cell(1,numel(vec));
    Hcoef = cell(1,numel(vec));
    for cnt = 1:numel(vec)
        cntb = vec(cnt);
        opt.ampErrBandWidth = cntb;
        [shape{cnt},comp{cnt},Hcoef{cnt}] = DeCom_MMD(sig,x,numGroup,ones(size(insAmp)),insFreq,insPhase,opt);
        
        pic = figure;
        plot((N1:N2)/Ntotal*480,comp{cnt}{1}(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
        pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-5,5]);
        xlabel('time'); ylabel('frequency');
        str = sprintf('./results/MMD_fig15_comp_bw_%d',cntb);
        saveas(pic,[str,'.fig']);
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
        
        
        pic = figure;
        plot((N1:N2)/Ntotal*480,sig(N1:N2)-comp{cnt}{1}(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
        pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-5,5]);
        str = sprintf('./results/MMD_fig15_comp_res_%d',cntb);
        saveas(pic,[str,'.fig']);
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    end
    save('./results/MMD_fig15.mat','-v7.3');
end

if (1)
    load ./results/MMD_fig15.mat;
    LL = length(shape{cnt}{1}.s0);
    pic = figure;plot(0:1/LL:(1-1/LL),shape{cnt}{1}.s0*Hcoef{cnt}{1}.a0);
    title('remove');axis square;
    set(gca, 'FontSize', 32);
    b=get(gca);
    set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
    tit = sprintf('./results/MMD_fig15_shape_%d.fig',0);
    saveas(pic,tit);
    str = sprintf('./results/MMD_fig15_shape_%d',0);
    print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
    
    LL = length(shape{cnt}{1}.scn{1});
    pic = figure;plot(0:1/LL:(1-1/LL),shape{cnt}{1}.scn{1}*Hcoef{cnt}{1}.a(1));
    title('remove');axis square;
    set(gca, 'FontSize', 32);
    b=get(gca);
    set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
    tit = sprintf('./results/MMD_fig15_shapeCn_%d.fig',1);
    saveas(pic,tit);
    str = sprintf('./results/MMD_fig15_shapeCn_%d',1);
    print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
    
    LL = length(shape{cnt}{1}.ssn{1});
    pic = figure;plot(0:1/LL:(1-1/LL),shape{cnt}{1}.ssn{1}*Hcoef{cnt}{1}.b(1));
    title('remove');axis square;
    set(gca, 'FontSize', 32);
    b=get(gca);
    set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
    tit = sprintf('./results/MMD_fig15_shapeSn_%d.fig',1);
    saveas(pic,tit);
    str = sprintf('./results/MMD_fig15_shapeSn_%d',1);
    print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
    
    LL = length(shape{cnt}{1}.scn{2});
    pic = figure;plot(0:1/LL:(1-1/LL),shape{cnt}{1}.scn{2}*Hcoef{cnt}{1}.a(2));
    title('remove');axis square;
    set(gca, 'FontSize', 32);
    b=get(gca);
    set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
    tit = sprintf('./results/MMD_fig15_shapeCn_%d.fig',-1);
    saveas(pic,tit);
    str = sprintf('./results/MMD_fig15_shapeCn_%d',-1);
    print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
    
    LL = length(shape{cnt}{1}.ssn{2});
    pic = figure;plot(0:1/LL:(1-1/LL),shape{cnt}{1}.ssn{2}*Hcoef{cnt}{1}.b(2));
    title('remove');axis square;
    set(gca, 'FontSize', 32);
    b=get(gca);
    set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
    tit = sprintf('./results/MMD_fig15_shapeSn_%d.fig',-1);
    saveas(pic,tit);
    str = sprintf('./results/MMD_fig15_shapeSn_%d',-1);
    print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
    
end




