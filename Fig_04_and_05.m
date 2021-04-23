% This code produces Figs. 4 and 5 from the paper "Decomposing non-stationary signals
% with time-varying wave-shape functions" by Marcelo A. Colominas and Hau-Tieng Wu.
%
% The running of this file might take several hours. Please, modify
% variable K (line 96) accordingly.
%
% Marcelo A. Colominas
% email: macolominas@conicet.gov.ar
% 23-APR-2021
addpath(genpath('DeCom-master')); %to compare with RDBR and MMD (developed by Haizhao Yang)

t = 0:0.001:1-0.001;%temporal variable
N = length(t);

phi1 = 2*(2*pi*3*t + 2*pi*3*t.^2);% phase function
phi2 = 2*(2*pi*5*t + 2*pi*3.5*t.^2 + 0.5*cos(2*pi*t));%phase function

%Defining the signals
f1 = 1.5*cos(phi1) + 0.25*cos(2.05*phi1) + 0.1*cos(3.05*phi1 + 0.01*phi1.^2) + 0.1*cos(4.05*phi1 + 0.01*phi1.^2); 
for i = 5:10
    f1 = f1 + 0.1*cos((i+0.05)*phi1+ 0.01*phi1.^2);
end;

f2 = cos(phi2);

for i = 2:10
    f2 = f2 + sqrt((1/i))*cos((i+0.01)*phi2);
end;

% Parameters for our methods ----------------------------------------
D1 = 10;
D2 = 10;
order_Phi = 3;
%--------------------------------------------------------------------

%-------Parameters for Estimating ridges and phases------------------

bt = 1:N;
redun = 4;
gamma = 0;
sigma = 0.25;
fmax = 0.05;
ft = 1/redun:1/redun:round(fmax*N);
%--------------------------------------------------------------------

%---Parameters for RDBR----------------------------------------------

N = length(t);

optR.maxiter = 10000;
optR.eps_error = 1e-12;
optR.show = 0;
optR.nknots = 20;
optR.knotremoval_factor= 1.01;
optR.order = 3;
optR.eps_diff = optR.eps_error;
    

ins_amplt = [ones(size(t));ones(size(t))];
    
fTrue = cell(1,2);
fTrue{1} = f1;
fTrue{2} = f2;
numGroup = 2;
fff = f1 + f2;
shapeTrue = cell(1,2);
%--------------------------------------------------------------------

% --- Parameters for MMD---------------------------------------------
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

% -------------------------------------------------------------------

index = 101:900;
SNR = [20 10 0 -5];

K = 50; %number of realizations
tiempo_LR = zeros(length(SNR),K);
tiempo_NLR = zeros(length(SNR),K);
tiempo_RDBR = zeros(length(SNR),K);
tiempo_MMD = zeros(length(SNR),K);

error1_LR = zeros(length(SNR),K);
error1_NLR = zeros(length(SNR),K);
error1_RDBR = zeros(length(SNR),K);
error1_MMD = zeros(length(SNR),K);

error2_LR = zeros(length(SNR),K);
error2_NLR = zeros(length(SNR),K);
error2_RDBR = zeros(length(SNR),K);
error2_MMD = zeros(length(SNR),K);

error1_tau = zeros(length(SNR),K);
error2_tau = zeros(length(SNR),K);



% Reference tau -----------------------------------------------------
[STFT,SST1,SST2,SST3,SST4,omega,omega2,omega3,omega4,tau2,tau3,phi22p,phi23p,phi24p] = sstn_test_modL_new(f1+f2,gamma,sigma,ft,bt,redun);
RTF = SST2;
tau_ref_1 = zeros(size(t));
tau_ref_2 = zeros(size(t));
jump = redun/2;
d = 4;
c1 = exridge(RTF,0,0,jump);
for i = 1:N
    tau_ref_1(i) = tau2(c1(i),i);
    RTF(c1(i)-d:c1(i)+d,i) = 0;
end;

c2 = exridge(RTF,0,0,jump);

for i = 1:N
    tau_ref_2(i) = tau2(c2(i),i);
end;
%--------------------------------------------------------------------

amp = 10.^(-SNR/20);
parfor j = 1:4
    tau_1 = zeros(size(t));
    tau_2 = zeros(size(t));
    aux = zeros(size(t));

for k = 1:K

f = f1+f2 + amp(j)*std(f1+f2)*randn(1,N); 

%----Estimating the ridges and phases and tau------------------------

[STFT,SST1,SST2,SST3,SST4,omega,omega2,omega3,omega4,tau2,tau3,phi22p,phi23p,phi24p] = sstn_test_modL_new(f,gamma,sigma,ft,bt,redun);

RTF = SST2;
c1 = exridge(RTF,0,0,jump);

for i = 1:N
    a = max(1,c1(i));
    a = min(c1(i),200);
    tau_1(i) = tau2(c1(i),i);
    a = max(1,c1(i)-d);
    b = min(200,c1(i)+d);
    aux(i) = sum(RTF(a:b,i));
    RTF(a:b,i) = 0;
end;
phi1_est = phase(aux);
A1_est = abs(aux);
c2 = exridge(RTF,0,0,jump);

for i = 1:N
    a = max(1,c2(i));
    a = min(c2(i),200);
    tau_2(i) = tau2(c2(i),i);
    a = max(1,c2(i)-d);
    b = min(200,c2(i)+d);
    aux(i) = sum(RTF(a:b,i));
end;
phi2_est = phase(aux);
A2_est = abs(aux);

%--------------------------------------------------------------------
error1_tau(j,k) = norm(tau_1(index)-tau_ref_1(index))/sqrt(800);
error2_tau(j,k) = norm(tau_2(index)-tau_ref_2(index))/sqrt(800);
%--------------------------------------------------------------------

ins_pre_phase = [phi1_est;phi2_est]*0.5/pi;
tic;
[shapeInLoop,comp_RDBR,errorRec,SL2Rec,iter,flag] = srcIterRegJC(f,N,numGroup,ins_amplt,ins_pre_phase,optR,fTrue,shapeTrue);
tiempo_RDBR(j,k) = toc;
error1_RDBR(j,k) = norm(comp_RDBR{1}(index)-f1(index))/sqrt(800);
error2_RDBR(j,k) = norm(comp_RDBR{2}(index)-f2(index))/sqrt(800);

tic;
warning('off')
[f_est,f1_est_NLR,f2_est_NLR,b_e] = two_wsf_nonlinear(f,ones(1,N),phi1_est,D1,ones(1,N),phi2_est,D2,order_Phi);
warning('on')
tiempo_NLR(j,k) = toc;
error1_NLR(j,k) = norm(f1_est_NLR(index)-f1(index))/sqrt(800);
error2_NLR(j,k) = norm(f2_est_NLR(index)-f2(index))/sqrt(800);

tic;
[f_est_lin,f1_est_LR,f2_est_LR,~,~,coef_est] = two_wsf_linear(f,ones(1,N),phi1_est,D1,ones(1,N),phi2_est,D2);
tiempo_LR(j,k) = toc;
error1_LR(j,k) = norm(f1_est_LR(index)-f1(index))/sqrt(800);
error2_LR(j,k) = norm(f2_est_LR(index)-f2(index))/sqrt(800);

tic;
inst_freq = zeros(2,N);
inst_freq(1,1:end-1) = N*diff(phi1_est)*0.5/pi; inst_freq(1,end) = inst_freq(1,end-1);
inst_freq(2,1:end-1) = N*diff(phi2_est)*0.5/pi; inst_freq(2,end) = inst_freq(2,end-1);
[shape,comp_MMD,Hcoef,flag,idx,iter,iterDR] = DeCom_MMD(f,t,2,[ones(size(t));ones(size(t))],inst_freq,[phi1_est;phi2_est]*0.5/pi,opt);
tiempo_MMD(j,k) = toc;
error1_MMD(j,k) = norm(comp_MMD{1}(index)-f1(index))/sqrt(800);
error2_MMD(j,k) = norm(comp_MMD{2}(index)-f2(index))/sqrt(800);
k
end;%k
end;%j

e1_MMD_media = zeros(1,4);
e2_MMD_media = zeros(1,4);
e1_MMD_std = zeros(1,4);
e2_MMD_std = zeros(1,4);

for i = 1:4
    e1_MMD_media(i) = mean(error1_MMD(i,error1_MMD(i,:)<2.5));
    e2_MMD_media(i) = mean(error2_MMD(i,error2_MMD(i,:)<2.5));
    e1_MMD_std(i) = std(error1_MMD(i,error1_MMD(i,:)<2.5));
    e2_MMD_std(i) = std(error2_MMD(i,error2_MMD(i,:)<2.5));
end;


% Graphics for Fig. 04-----------------------------------------------

figure;subplot(221)
errorbar(([0 1 2 2.5]),mean(error1_RDBR(:,:),2), std(error1_RDBR(:,:)'),'s-'); hold on
errorbar(([0 1 2 2.5]),mean(error1_NLR(:,:),2), std(error1_NLR(:,:)'),'*-'); 
errorbar(([0 1 2 2.5]),mean(error1_LR(:,:),2),std(error1_LR(:,:)'),'o-'); 
errorbar(([0 1 2 2.5]),e1_MMD_media,e1_MMD_std,'^-');
legend('RDBR','SAMD','LR','MMD','Location','Northwest')
ylabel('RMSE')
xlabel('SNR in (dB)')
xlim([-0.1 2.6])
title('Error for $s_1(t)$','interpreter','latex')
set(gca,'xtick',[0 1 2 2.5])
set(gca,'xticklabel',SNR)

subplot(222)
errorbar(([0 1 2 2.5]),mean(error2_RDBR(:,:),2),std(error2_RDBR(:,:)'),'s-'); hold on
errorbar(([0 1 2 2.5]),mean(error2_NLR(:,:),2),std(error2_NLR(:,:)'),'*-'); 
errorbar(([0 1 2 2.5]),mean(error2_LR(:,:),2),std(error2_LR(:,:)'),'s-'); 
errorbar(([0 1 2 2.5]),e2_MMD_media,e2_MMD_std,'^-');
ylabel('RMSE')
legend('RDBR','SAMD','LR','MMD','Location','Northwest')
xlabel('SNR in (dB)')
xlim([-0.1 2.6])
title('Error for $s_2(t)$','interpreter','latex')
set(gca,'xtick',[0 1 2 2.5])
set(gca,'xticklabel',SNR)


%Graphics for Fig. 05 -----------------------------------------------
figure; subplot(221)
scatter(error1_NLR(1,:),error1_tau(1,:));
hold on;
scatter(error1_NLR(2,:),error1_tau(2,:));
scatter(error1_NLR(3,:),error1_tau(3,:));
scatter(error1_NLR(4,:),error1_tau(4,:));
box on
legend('20 dB','10 dB','0 dB', '-5 dB')
xlabel('RMSE of $s_1(t)$ (SAMD)','interpreter','latex')
ylabel('RMSE of group delay','interpreter','latex')
subplot(222)
scatter(error2_NLR(1,:),error2_tau(1,:));
hold on;
scatter(error2_NLR(2,:),error2_tau(2,:));
scatter(error2_NLR(3,:),error2_tau(3,:));
scatter(error2_NLR(4,:),error2_tau(4,:));
box on
xlim([0.4 1.5])
xlabel('RMSE of $s_2(t)$ (SAMD)','interpreter','latex')

%--------------------------------------------------------------------

