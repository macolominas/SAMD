% This code generate Figure 6 in the paper "Recursive Diffeomorphism-Based
% Regression for Shape Functions".
%
% By Haizhao Yang

clear all;
close all;

if(1)
    %stepSize = 1/nbins;
    nbins = 50;
    
    %set up data
    N = 8192/2;
    x = [0:N-1]/N;
    amp = 0.01;
    F1 = 60;
    F2 = 90;
    xx = x + amp*sin(2*pi*x);
    yy = x + amp*cos(2*pi*x);
    
    ins_pre_phase_1 = F1*xx;
    ins_pre_phase_2 = F2*yy;
    
    X1 = mod(ins_pre_phase_1, 1);
    X2 = mod(ins_pre_phase_2, 1);
   
    X = rand(1,N);
    pic = figure;
    h = histogram(X,nbins); saveas(pic,'./results/RDBR_fig6_1.fig'); axis square;
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDBR_fig6_1';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    
    pic = figure;
    h1 = histogram(X1,nbins); saveas(pic,'./results/RDBR_fig6_2.fig'); axis square;
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDBR_fig6_2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    
    pic = figure;
    h2 = histogram(X2,nbins); saveas(pic,'./results/RDBR_fig6_3.fig'); axis square;
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDBR_fig6_3';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    
    
    %set up data
    N = 8192*2;
    x = [0:N-1]/N;
    amp = 0.01;
    F1 = 60;
    F2 = 90;
    xx = x + amp*sin(2*pi*x);
    yy = x + amp*cos(2*pi*x);
    
    ins_pre_phase_1 = F1*xx;
    ins_pre_phase_2 = F2*yy;
    
    X1 = mod(ins_pre_phase_1, 1);
    X2 = mod(ins_pre_phase_2, 1);
    
    %stepSize = 1/nbins;
    X = rand(1,N);
    pic = figure;
    h = histogram(X,nbins); saveas(pic,'./results/RDBR_fig6_4.fig'); axis square;
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDBR_fig6_4';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    
    pic = figure;
    h1 = histogram(X1,nbins); saveas(pic,'./results/RDBR_fig6_5.fig'); axis square;
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDBR_fig6_5';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    
    pic = figure;
    h2 = histogram(X2,nbins); saveas(pic,'./results/RDBR_fig6_6.fig'); axis square;
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDBR_fig6_6';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    
    
    %set up data
    N = 8192*8;
    x = [0:N-1]/N;
    amp = 0.01;
    F1 = 60;
    F2 = 90;
    xx = x + amp*sin(2*pi*x);
    yy = x + amp*cos(2*pi*x);
    
    ins_pre_phase_1 = F1*xx;
    ins_pre_phase_2 = F2*yy;
    
    X1 = mod(ins_pre_phase_1, 1);
    X2 = mod(ins_pre_phase_2, 1);
    
    %stepSize = 1/nbins;
    X = rand(1,N);
    pic = figure;
    h = histogram(X,nbins); saveas(pic,'./results/RDBR_fig6_7.fig'); axis square;
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDBR_fig6_7';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    
    pic = figure;
    h1 = histogram(X1,nbins); saveas(pic,'./results/RDBR_fig6_8.fig'); axis square;
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDBR_fig6_8';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    
    pic = figure;
    h2 = histogram(X2,nbins); saveas(pic,'./results/RDBR_fig6_9.fig'); axis square;
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDBR_fig6_9';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    
    
end