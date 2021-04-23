% This code generate Figure 1 in the paper "Recursive Diffeomorphism-Based
% Regression for Shape Functions".
%
% By Haizhao Yang

clear all;
close all;

load 'ECGshape2.mat';
N = length(y);
x = 0:1/N:(1-1/N);
z=fftshift(fft(ifftshift(y)))/length(y);
%first 10
z2=zeros(size(z));
t = 500-10:500+10;
z2(t)=z(t);
y2=real(fftshift(ifft(ifftshift(z2)))*length(y));
pic = figure;hold on;
plot(x,y,'b-','LineWidth',4);
plot(x,y2,'r.','LineWidth',4);
axis square;
xlabel('t');%ylabel('s(t)');
set(gca, 'FontSize', 16);
b=get(gca);
set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
tit = sprintf('./results/RDBR_fig1_%d.fig',1);
saveas(pic,tit);
str = sprintf('./results/RDBR_fig1_%d',1);
print(gcf, '-depsc', str);

pic = figure;plot(-100:100,abs(z(413:613)),'LineWidth',4);
axis square;
set(gca, 'FontSize', 16);
b=get(gca);
set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
tit = sprintf('./results/RDBR_fig1_%d.fig',2);
saveas(pic,tit);
str = sprintf('./results/RDBR_fig1_%d',2);
print(gcf, '-depsc', str);