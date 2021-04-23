% This code provides an example in which the signal consists of two MIMFs.
% We use MMD and GMD methods for the decomposition to see the difference of
% the performance. We will see that GMD fails to capture some information
% of MIMFs and hence the residual signal still contain significant
% oscillation patterns.
%
% Note: Please run MMD_Fig1_3.m first.
% 
% Reference:
% "Multiresolution Mode Decomposition (MMD) for Adaptive Time Series Analysis"
% by Haizhao Yang.

if (1)
    close all;
    clear all;
    
    %% generate signal
    load ./results/MMD_fig1.mat;
    sig = Sig; % Clean MMD model with two components from PPG signals.
    
    x = [0:N-1]/N;
    
    numGroup = 2;
    opt.eps = 1e-3;
    opt.res = 0.25;
    opt.freq_range = [0 N/32];
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
    
end

% perform MMD
if (1)
    opt.maxiter = 200;
    opt.eps_error = 1e-6;
    opt.show = 0;
    opt.iterStyle = 'GS';
    opt.shapeMethod = 2;
    opt.eps_diff = 1e-6;
    opt.ampErrBandWidth = 40;
    opt.numSweep = 10;
    
    switch opt.shapeMethod
        case 1
            opt.para.Ls=1000;
            opt.para.bandWidth = 10;
            opt.para.diffeoMethod = 'nufft';
        case 2
            opt.para.nknots = 10;
            opt.para.knotremoval_factor= 1.0001;
            opt.para.order = 3;
            opt.para.Ls = 1000;
    end
    
    % test example: two components
    [shape,comp] = DeCom_MMD(sig,x,numGroup,ones(size(insAmp)),insFreq,insPhase,opt);
    
    %% Modified to regenerate the signal to test if it is actually accurate enough
    Sig = comp{1} + comp{2};
    save('./results/MMD_fig21_MMD.mat','-v7.3');
end

if (1)
    load ./results/MMD_fig21_MMD.mat;
    close all;
    N1 = 1; N2 = N;
    Ntotal = length(signal.pleth.y);
    pic = figure;
    plot((N1:N2)/Ntotal*480,comp{1}(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    title('First component by MMD');
    saveas(pic,'./results/MMD_fig21_comp1_MMD.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig21_comp1_MMD';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,comp{2}(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    title('Second component by MMD');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/MMD_fig21_comp2_MMD.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig21_comp2_MMD';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,sig(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    title('Original signal');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/MMD_fig21_org_MMD.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig21_org_MMD';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,sig(N1:N2)-Sig(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    title('Error of MMD');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/MMD_fig21_res_MMD.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig21_res_MMD';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end

% perform GMD
if (1)
    opt.maxiter = 200;
    opt.eps_error = 1e-6;
    opt.show = 0;
    opt.iterStyle = 'GS';
    opt.shapeMethod = 2;
    opt.eps_diff = 1e-6;
    opt.ampErrBandWidth = 0;
    opt.numSweep = 1;
    
    switch opt.shapeMethod
        case 1
            opt.para.Ls=1000;
            opt.para.bandWidth = 10;
            opt.para.diffeoMethod = 'nufft';
        case 2
            opt.para.nknots = 10;
            opt.para.knotremoval_factor= 1.0001;
            opt.para.order = 3;
            opt.para.Ls = 1000;
    end
    
    % test example: two components
    [shape1,comp1] = DeCom_MMD(sig,x,numGroup,insAmp,insFreq,insPhase,opt);
    
    %% Modified to regenerate the signal to test if it is actually accurate enough
    Sig1 = comp1{1} + comp1{2};
    save('./results/MMD_fig21_GMD.mat','-v7.3');
end

if (1) 
    load ./results/MMD_fig21_GMD.mat;
    close all;
    N1 = 1; N2 = N;
    Ntotal = length(signal.pleth.y);
    pic = figure;
    plot((N1:N2)/Ntotal*480,comp1{1}(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    title('First component by GMD');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/MMD_fig21_comp1_p2_GMD.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig21_comp1_p2_GMD';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,comp1{2}(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    title('Second component by GMD');
    saveas(pic,'./results/MMD_fig21_comp2_p2_GMD.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig21_comp2_p2_GMD';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,sig(N1:N2)-Sig1(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    title('Error of GMD');
    saveas(pic,'./results/MMD_fig21_res_p2_GMD.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig21_res_p2_GMD';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end
