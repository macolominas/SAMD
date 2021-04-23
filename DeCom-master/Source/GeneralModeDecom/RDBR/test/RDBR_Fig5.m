% This code generate Figure 5 in the paper "Recursive Diffeomorphism-Based
% Regression for Shape Functions".
%
% By Haizhao Yang

clear all;
close all;

if(1)
    %set up data
    N = 8192/4;
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
    
    insPhase = F1*xx;
    insAmplitude = am;
    
    Y = gen_shape2(x,2);
    pic = figure; plot(x,Y,'.','LineWidth',4);axis square;saveas(pic,'./results/RDBR_fig4_1.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDBR_fig4_1';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    phi2 = mod(insPhase, 1);
    X = phi2(:);   Y = fff./insAmplitude; Y = Y(:);
    pic = figure; plot(X,Y,'.','LineWidth',4);axis square;saveas(pic,'./results/RDBR_fig4_2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDBR_fig4_2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    options = struct('animation', 0, 'knotremoval_factor',1.0001);
    [pp] = BSFK(X',Y,3,20,[],options);
    x = 0:1/1000:1;
    shapeEst = ppval(pp,x);
    shapeTrue = gen_shape2(x,2);
    pic = figure; 
    h(1) = plot(x,shapeTrue,'r','LineWidth',2); hold on;h(2) = plot(x,shapeEst,'b','LineWidth',2); 
    legend(h,'true','estimated','northeast');
    axis square;saveas(pic,'./results/RDBR_fig3_4.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDBR_fig4_3';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end