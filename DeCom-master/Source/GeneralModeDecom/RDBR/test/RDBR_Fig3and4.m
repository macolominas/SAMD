% This code generate Figure 3 and 4 in the paper "Recursive Diffeomorphism-Based
% Regression for Shape Functions".
%
% By Haizhao Yang

clear all;
close all;

if(1)
    %set up data
    N = 8192*2;
    x = [0:N-1]/N;
    fff = zeros(1,N);
    amp = 0.01;
    F1 = 60;
    F2 = 90;
    xx = x + amp*sin(2*pi*x);
    f1 = zeros(1,N);
    am = 1+0.05*sin(4*pi*x);
    f1 = am.*gen_shape2(F1*xx,2);
    yy = x + amp*cos(2*pi*x);
    %yy = x + amp*sin(4*pi*x);
    f2 = zeros(1,N);
    bm = 1+0.1*sin(2*pi*x);
    f2 = bm.*gen_shape2(F2*yy,3);
    
    NM = 0;% is fine
    ns = NM*randn(1,N);
    fff = f1 + f2 + ns;
    SNR = min(10*log10(var(f1)/NM^2),10*log10(var(f2)/NM^2))
    
    
    z = gen_shape2(x,2);
    zz = abs(fftshift(fft(ifftshift(z)))/length(z));
    pic = figure; plot(x,z,'LineWidth',4);axis square;saveas(pic,'./results/RDBR_fig2_1.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDBR_fig2_1';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    
    pic = figure; plot(-50:50,zz(N/2-49:N/2+51),'LineWidth',4);axis square;saveas(pic,'./results/RDBR_fig2_2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDBR_fig2_2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    z = gen_shape2(x,3);
    zz = abs(fftshift(fft(ifftshift(z)))/length(z));
    pic = figure; plot(x,z,'LineWidth',4);axis square;saveas(pic,'./results/RDBR_fig2_3.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDBR_fig2_3';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure; plot(-50:50,zz(N/2-49:N/2+51),'LineWidth',4);axis square;saveas(pic,'./results/RDBR_fig2_4.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDBR_fig2_4';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end

if (1)
    %synchrosqueezed wave packet transform
    eps = 1e-3;
    res = 0.5;
    freq_range = [0 N/2];
    NG = N/2;
    dt = 1/NG;
    red = 1;
    t_sc = 0.7;
    rad = 1;
    [T_f coef kk] = ss_wp1_fwd(fff,1,1,1,NG,x,freq_range(2),freq_range(1),rad,1,t_sc,red,eps,res);
    loc = find(T_f<0.05*NM^2);
    T_f(loc) = 0;
    pic = figure;imagesc([0 1],[1 150],(T_f(1:300,:)));colormap (1-gray);axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    tit = sprintf('./results/RDBR_fig3_%d.fig',1);
    saveas(pic,tit);
    str = sprintf('./results/RDBR_fig3_%d',1);
    print(gcf, '-depsc', str);
    
    %synchrosqueezed wave packet transform
    eps = 1e-3;
    res = 2;
    freq_range = [0 N/2];
    NG = N/2;
    dt = 1/NG;
    red = 1;
    t_sc = 0.7;
    rad = 1;
    [T_f coef kk] = ss_wp1_fwd(fff,1,1,1,NG,x,freq_range(2),freq_range(1),rad,1,t_sc,red,eps,res);
    loc = find(T_f<0.05*NM^2);
    T_f(loc) = 0;
    pic = figure;imagesc([0 1],[1 N/8],log10(T_f(1:N/16,:)));colormap (1-gray);axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    tit = sprintf('./results/RDBR_fig3_%d.fig',2);
    saveas(pic,tit);
    str = sprintf('./results/RDBR_fig3_%d',2);
    print(gcf, '-depsc', str);
end
