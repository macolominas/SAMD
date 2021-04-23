% This code produces Fig. 19 of the Supplemental Material from the paper "Decomposing non-stationary signals
% with time-varying wave-shape functions" by Marcelo A. Colominas and Hau-Tieng Wu.
%
%
% Marcelo A. Colominas
% email: macolominas@conicet.gov.ar
% 23-APR-2021
addpath(genpath('DeCom-master')); %to compare with RDBR and MMD (developed by Haizhao Yang)

x = load('f1o07m.mat');
x = x.val;

x = x(:)';
x = x - mean(x);
fs = 250;
N = length(x);
t = 0:1/fs:N/fs-1/fs;


% Parameters for STFT
bt = 1:N;
redun = 1;
gamma = 0;
sigma = 0.15;
fmax = 0.05;
ft = 1/redun:1/redun:round(fmax*N);

[STFT,SST1,SST2,SST3,SST4,omega,omega2,omega3,omega4,tau2,tau3,phi22p,phi23p,phi24p] = sstn_test_modL(x,gamma,sigma,ft,bt);

RTF = SST2;
d = 2;
c1 = exridge(RTF(1:16,:),0,0,2);
aux = zeros(size(t));
for j = 1:N
    aux(j) = sum(RTF(c1(j)-d:c1(j)+d,j));
    RTF(c1(j)-d:c1(j)+d,j) = 0;
end;
A1_est = abs(aux); 
phi1_est = phase(aux);

D1 = 40; 

tic;
C = [cosenos(A1_est,phi1_est,D1)];

C1_aux = cosenos(ones(1,1000),-pi:2*pi/1000:pi-2*pi/1000,D1);

coef_est = (x*C')*inv(C*C');
f_est_LR = coef_est*C;
f1_est_LR = coef_est(1:2*D1)*C(1:2*D1,:);
tiempo_LR = toc

WSF1 = coef_est(1:2*D1)*C1_aux;

order_Phi = 1;
warning('off')
tic;
[f1_est_NLR,b_e] = one_wsf_nonlinear(x,A1_est,phi1_est,D1,order_Phi);
tiempo_NLR = toc
warning('on')

N = length(t);

opt.maxiter = 10000;
opt.eps_error = 1e-12;
opt.show = 0;
opt.nknots = 20;
opt.knotremoval_factor= 1.01;
opt.order = 3;
opt.eps_diff = opt.eps_error;
    
ins_pre_phase = [phi1_est]*0.5/pi;
ins_amplt = [A1_est];
    
fTrue = cell(1,2);
numGroup = 1;

shapeTrue = cell(1,2);
%         shapeTrue{1} = @(x) sh1(x);
%         shapeTrue{2} = @(x) sh2(x);
tic;
[shapeInLoop,comp,errorRec,SL2Rec,iter,flag] = srcIterRegJC(x,N,numGroup,ins_amplt,ins_pre_phase,opt,fTrue,shapeTrue);
tiempo_RDBR = toc

opt.maxiter = 300;
opt.eps_error = 1e-6;
opt.show = 0;
opt.iterStyle = 'GS';
opt.shapeMethod = 2;
opt.eps_diff = 1e-6;
opt.ampErrBandWidth = 20;
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
tic;
inst_freq = zeros(2,N);
inst_freq(1,1:end-1) = N*diff(phi1_est)*0.5/pi; inst_freq(1,end) = inst_freq(1,end-1);
[shape,component,Hcoef,flag,idx,iter,iterDR] = DeCom_MMD(x,t,1,[A1_est],inst_freq,[phi1_est]*0.5/pi,opt);
tiempo_MMD = toc

%----Graphics -------------------------------------------------------

%--------------------------------------------------------------------
ind = 201:2400;
N = length(f1_est_NLR(ind));
t_m = linspace(min(phi1_est(ind)),max(phi1_est(ind)),N);
T = length(WSF1)/5;
P = ((t_m(end)-t_m(1)+1)/(2*pi));
P = floor(P);
m = interp1(phi1_est(ind),1:N,t_m,'spline');

y = interp1(1:N,f1_est_NLR(ind)./A1_est(ind),m,'spline');
waves_SAMD = zeros(P,T); 
for i = 1:P
    waves_SAMD(i,:) = interp1(0.5*t_m/pi,y,0.5*t_m(1)/pi+i-1+[0:1/(T):1-1/(T)],'spline');
end;

y = interp1(1:N,comp{1}(ind)./A1_est(ind),m,'spline');
waves_RDBR = zeros(P,T); 
for i = 1:P
    waves_RDBR(i,:) = interp1(0.5*t_m/pi,y,0.5*t_m(1)/pi+i-1+[0:1/(T):1-1/(T)],'spline');
end;

y = interp1(1:N,component{1}(ind)./A1_est(ind),m,'spline');
waves_MMD = zeros(P,T); 
for i = 1:P
    waves_MMD(i,:) = interp1(0.5*t_m/pi,y,0.5*t_m(1)/pi+i-1+[0:1/(T):1-1/(T)],'spline');
end;



figure;
subplot(5,5,1:4);
plot(t(ind)-t(ind(1)),x(ind),'k')
text(0.05,6500,'ECG recording f1o07m')
ylim([-1800 7200])

subplot(5,5,6:9);
plot(t(ind)-t(ind(1)),f1_est_NLR(ind),'k');
text(0.05,4500,'SAMD')
subplot(5,5,10)
plot(linspace(-0.5,0.5,T),waves_SAMD')
ylim([-1.8 7.1])

subplot(5,5,11:14);
plot(t(ind)-t(ind(1)),f1_est_LR(ind),'k');
text(0.05,3500,'LR')
subplot(5,5,15);
plot(linspace(-0.5,0.5,length(WSF1)),WSF1)
ylim([-1.8 7.1])

subplot(5,5,16:19)
plot(t(ind)-t(ind(1)),comp{1}(ind),'k');
text(0.05,2800,'RDBR')
ylim([-815 3200])
subplot(5,5,20)
plot(linspace(-0.5,0.5,T),waves_RDBR')
ylim([-1.8 7.1])

subplot(5,5,21:24);
plot(t(ind)-t(ind(1)),component{1}(ind),'k');
text(0.05,5500,'MMD')
ylim([-1600 6000])
xlabel('t(s)','interpreter','latex')
subplot(5,5,25)
plot(linspace(-0.5,0.5,T),waves_MMD')
ylim([-1.8 7.1])

