% This code produces Fig. 7 from the paper "Decomposing non-stationary signals
% with time-varying wave-shape functions" by Marcelo A. Colominas and Hau-Tieng Wu.
%
%
% Marcelo A. Colominas
% email: macolominas@conicet.gov.ar
% 23-APR-2021
addpath(genpath('DeCom-master')); %to compare with RDBR and MMD (developed by Haizhao Yang)

x = load('2869156_Resp(32).txt');
ecg = load('2869156_II(256).txt');
pleth = load('2869156_Pleth(64).txt');

x = x(:)';
fs = 32;
x = x(1:fs*60*10);
N = length(x);
t = 0:1/fs:N/fs-1/fs;

fs_ecg = 256;
ecg = ecg(1:fs_ecg*60*10);
t_ecg = 0:1/fs_ecg:length(ecg)/fs_ecg-1/fs_ecg;

fs_pleth = 64;
pleth = pleth(1:fs_pleth*60*10)';
t_pleth = 0:1/fs_pleth:length(pleth)/fs_pleth-1/fs_pleth;


% Parameters for STFT and SST ---------------------------------------
bt = 1:N;
redun = 1;
gamma = 0;
sigma = 0.045;
fmax = 0.05;
ft = 1/redun:1/redun:round(fmax*N);

[STFT,SST1,SST2,SST3,SST4,omega,omega2,omega3,omega4,tau2,tau3,phi22p,phi23p,phi24p] = sstn_test_modL(x,gamma,sigma,ft,bt);

% Estimating ridges, amplitudes and phases --------------------------
SST2(1:20,:) = 0; STFT(1:20,:) = 0;
RTF = SST2;
d = 2;
c1 = exridge(RTF,0,0,2);
aux = zeros(size(t));
for j = 1:N
    aux(j) = sum(RTF(c1(j)-d:c1(j)+d,j));
    RTF(c1(j)-d:c1(j)+d,j) = 0;
end;
A1_est = abs(aux); 
phi1_est = phase(aux);

c2 = exridge(RTF(201:350,:),0,0,2)+200;
aux = zeros(size(t));
for j = 1:N
    aux(j) = sum(RTF(c2(j)-d:c2(j)+d,j));
end;
A2_est = abs(aux); 
phi2_est = phase(aux);

D1 = 2; 
D2 = 5; 

tic;
C = [cosenos(A1_est,phi1_est,D1);cosenos(A2_est,phi2_est,D2)];

C1_aux = cosenos(ones(1,1000),-pi:2*pi/1000:pi-2*pi/1000,D1);
C2_aux = cosenos(ones(1,1000),-pi:2*pi/1000:pi-2*pi/1000,D2);

coef_est = (x*C')*inv(C*C');
f_est_LR = coef_est*C;
f1_est_LR = coef_est(1:2*D1)*C(1:2*D1,:);
f2_est_LR = coef_est(2*D1+1:end)*C(2*D1+1:end,:);
tiempo_LR = toc

WSF1 = coef_est(1:2*D1)*C1_aux;
WSF2 = coef_est(2*D1+1:end)*C2_aux;

order_Phi = 2;
warning('off')
tic;
[f_est_NLR,f1_est_NLR,f2_est_NLR,b_e] = two_wsf_nonlinear(x,A1_est,phi1_est,D1,A2_est,phi2_est,D2,order_Phi);
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
    
ins_pre_phase = [phi1_est;phi2_est]*0.5/pi;
ins_amplt = [A1_est;A2_est];
    
fTrue = cell(1,2);
numGroup = 2;

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
inst_freq(2,1:end-1) = N*diff(phi2_est)*0.5/pi; inst_freq(2,end) = inst_freq(2,end-1);
[shape,component,Hcoef,flag,idx,iter,iterDR] = DeCom_MMD(x,t,2,[A1_est;A2_est],inst_freq,[phi1_est;phi2_est]*0.5/pi,opt);
tiempo_MMD = toc

%----Graphics -------------------------------------------------------

ind = 180*fs:210*fs;
t_m = linspace(min(phi2_est(ind)),max(phi2_est(ind)),N);
N = length(f2_est_NLR(ind));
T = length(WSF2);
P = ((t_m(end)-t_m(1)+1)/(2*pi));
P = floor(P);
m = interp1(phi2_est(ind),1:N,t_m,'spline');

y = interp1(1:N,f2_est_NLR(ind)./A2_est(ind),m,'spline');
waves_SAMD = zeros(P,T); 
for i = 1:P
    waves_SAMD(i,:) = interp1(0.5*t_m/pi,y,0.5*t_m(1)/pi+i-1+[0:1/(T):1-1/(T)],'spline');
end;

y = interp1(1:N,comp{2}(ind)./A2_est(ind),m,'spline');
waves_RDBR = zeros(P,T);
for i = 1:P
    waves_RDBR(i,:) = interp1(0.5*t_m/pi,y,0.5*t_m(1)/pi+i-1+[0:1/(T):1-1/(T)],'spline');
end;

y = interp1(1:N,component{2}(ind)./A2_est(ind),m,'spline');
waves_MMD = zeros(P,T);
for i = 1:P
    waves_MMD(i,:) = interp1(0.5*t_m/pi,y,0.5*t_m(1)/pi+i-1+[0:1/(T):1-1/(T)],'spline');
end;

figure;subplot(5,9,1:4);
ind = 180*fs:210*fs;
ind_ecg = 180*fs_ecg:210*fs_ecg;
plot(t(ind),x(ind),'k')
text(181,2800,'IP recording')
%ylim([-10 15])
xlim([180 210])

subplot(5,9,10:13)
plot(t(ind),f1_est_NLR(ind),'k')
text(181,900,'SAMD')
% ylim([-6.5 6.5])
xlim([180 210])

subplot(5,9,14:17)
plot(t_ecg(ind_ecg),ecg(ind_ecg)-8200,'r-.');hold on; plot(t(ind),f2_est_NLR(ind),'k'); 
text(181,240,'SAMD')
%ylim([-5 9])
xlim([180 210])

subplot(5,9,18)
for i=[1 P]
plot(linspace(-0.5,0.5,length(WSF2)),waves_SAMD(i,:)); hold on
end;
ylim([0.02 0.14])

% set(gca,'yticklabel',[])

subplot(5,9,19:22)
plot(t(ind),f1_est_LR(ind),'k'); 
text(181,900,'LR')
% ylim([-6.5 6.5])
xlim([180 210])

subplot(5,9,23:26)
plot(t_ecg(ind_ecg),ecg(ind_ecg)-8200,'r-.');hold on; plot(t(ind),f2_est_LR(ind),'k');
text(181,240,'LR')
% ylim([-5 9])
xlim([180 210])

subplot(5,9,27);
plot(linspace(-0.5,0.5,length(WSF2)),WSF2)
ylim([0.02 0.14])

% set(gca,'yticklabel',[])

subplot(5,9,28:31)
plot(t(ind),comp{1}(ind),'k');
text(181,1500,'RDBR')
% ylim([-6.5 6.5])
xlim([180 210])

subplot(5,9,32:35)
plot(t_ecg(ind_ecg),ecg(ind_ecg)-8200,'r-.');hold on; plot(t(ind),comp{2}(ind),'k');
text(181,300,'RDBR')
% ylim([-5 9])
xlim([180 210])

subplot(5,9,36)
plot(linspace(-0.5,0.5,length(WSF2)),waves_RDBR')
% set(gca,'yticklabel',[])

subplot(5,9,37:40)
plot(t(ind),component{1}(ind),'k');
text(181,1000,'MMD')
% ylim([-11 11])
xlim([180 210])
xlabel('$t$ (s)','interpreter','latex')

subplot(5,9,41:44)
plot(t_ecg(ind_ecg),ecg(ind_ecg)-8200,'r-.');hold on; plot(t(ind),component{2}(ind),'k');
text(181,700,'MMD')
% ylim([-5 9])
xlim([180 210])
xlabel('$t$ (s)','interpreter','latex')

subplot(5,9,45)
for i=1:P
plot(linspace(-0.5,0.5,length(WSF2)),waves_MMD(i,:)); hold on
end;
% set(gca,'yticklabel',[])

