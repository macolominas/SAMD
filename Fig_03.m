% This code produces Fig. 3 from the paper "Decomposing non-stationary signals
% with time-varying wave-shape functions" 
% by Marcelo A. Colominas and Hau-Tieng Wu.
%
% The running of this file might take several hours. Please, modify
% variable K (line 95) accordingly.
%
% Marcelo A. Colominas
% email: macolominas@conicet.gov.ar
% 23-APR-2021
addpath(genpath('DeCom-master')); %to compare with RDBR and MMD (developed by Haizhao Yang)
t = 0:0.001:1-0.001;% temporal variable
N = length(t);

phi1 = 2*(2*pi*3*t + 2*pi*3*t.^2);% phase function
phi2 = 2*(2*pi*5*t + 2*pi*3.5*t.^2 + 0.5*cos(2*pi*t));%phase function

% Defining the signals
f1 = 1.5*cos(phi1) + 0.25*cos(2.05*phi1) + 0.1*cos(3.05*phi1 + 0.01*phi1.^2) + 0.1*cos(4.05*phi1 + 0.01*phi1.^2); 
for i = 5:10
    f1 = f1 + 0.1*cos((i+0.05)*phi1+ 0.01*phi1.^2);
end;

f2 = cos(phi2);

for i = 2:10
    f2 = f2 + sqrt((1/i))*cos((i+0.01)*phi2);
end;

% Parameters for our methods ---------------------------------------
D1 = 10;
D2 = 10;
order_Phi = 3;
%-------------------------------------------------------------------

%-------Parameters for Estimating ridges and phases-----------------
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
    
ins_pre_phase = [phi1;phi2]*0.5/pi;
ins_amplt = [ones(size(t));ones(size(t))];
    
fTrue = cell(1,2);
fTrue{1} = f1;
fTrue{2} = f2;
numGroup = 2;
fff = f1 + f2;
shapeTrue = cell(1,2);
%         shapeTrue{1} = @(x) sh1(x);
%         shapeTrue{2} = @(x) sh2(x);
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

K = 100; %number of realizations
tiempo_LR = zeros(1,K);
tiempo_NLR = zeros(1,K);
tiempo_RDBR = zeros(1,K);
tiempo_MMD = zeros(1,K);

f1_est_LR = zeros(K,N);
f2_est_LR = zeros(K,N);
f1_est_NLR = zeros(K,N);
f2_est_NLR = zeros(K,N);

comp_RDBR = cell(1,K);
comp_MMD = cell(1,K);

f1_est_RDBR = zeros(K,N);
f2_est_RDBR = zeros(K,N);
f1_est_MMD = zeros(K,N);
f2_est_MMD = zeros(K,N);


parfor k = 1:K

f = f1+f2 + 0.3162*std(f1+f2)*randn(1,N);

%----Estimating the ridges and phases--------------------------------

[STFT,SST1,SST2,SST3,SST4,omega,omega2,omega3,omega4,tau2,tau3,phi22p,phi23p,phi24p] = sstn_test_modL_new(f,gamma,sigma,ft,bt,redun);


RTF = SST2;
jump = redun/2;
d = 4;
c1 = exridge(RTF,0,0,jump);
aux = zeros(1,N);

for i = 1:N
    a = max(1,c1(i)-d);
    b = min(200,c1(i)+d);
    aux(i) = sum(RTF(a:b,i));
    RTF(a:b,i) = 0;
end;
RTF(1:10,:) = 0;
phi1_est = phase(aux);
A1_est = abs(aux);
c2 = exridge(RTF,0,0,jump);

for i = 1:N
    a = max(1,c2(i)-d);
    b = min(200,c2(i)+d);
    aux(i) = sum(RTF(a:b,i));
end;
phi2_est = phase(aux);
A2_est = abs(aux);

%--------------------------------------------------------------------


ins_pre_phase = [phi1_est;phi2_est]*0.5/pi;
tic;
[shapeInLoop,comp_RDBR{k},errorRec,SL2Rec,iter,flag] = srcIterRegJC(f,N,numGroup,ins_amplt,ins_pre_phase,optR,fTrue,shapeTrue);
tiempo_RDBR(k) = toc;
f1_est_RDBR(k,:) = comp_RDBR{k}{1};
f2_est_RDBR(k,:) = comp_RDBR{k}{2};

tic;
warning('off')
[f_est,f1_est_NLR(k,:),f2_est_NLR(k,:),b_e] = two_wsf_nonlinear(f,ones(1,N),phi1_est,D1,ones(1,N),phi2_est,D2,order_Phi);
warning('on')
tiempo_NLR(k) = toc;

tic;
[f_est_lin,f1_est_LR(k,:),f2_est_LR(k,:),~,~,coef_est] = two_wsf_linear(f,ones(1,N),phi1_est,D1,ones(1,N),phi2_est,D2);
tiempo_LR(k) = toc;

tic;
inst_freq = zeros(2,N);
inst_freq(1,1:end-1) = N*diff(phi1_est)*0.5/pi; inst_freq(1,end) = inst_freq(1,end-1);
inst_freq(2,1:end-1) = N*diff(phi2_est)*0.5/pi; inst_freq(2,end) = inst_freq(2,end-1);
[shape,comp_MMD{k},Hcoef,flag,idx,iter,iterDR] = DeCom_MMD(f,t,2,[ones(size(t));ones(size(t))],inst_freq,[phi1_est;phi2_est]*0.5/pi,opt);
tiempo_MMD(k) = toc;
f1_est_MMD(k,:) = comp_MMD{k}{1};
f2_est_MMD(k,:) = comp_MMD{k}{2};
k
end;

disp('mean time NLR: ')
mean(tiempo_NLR)

disp('mean time LR: ')
mean(tiempo_LR)

disp('mean time RDBR: ')
mean(tiempo_RDBR)

disp('mean time MMD: ')
mean(tiempo_MMD)

index = 101:900;

disp('mean and std f1 RMSE NLR:')
mean(sqrt(sum( (f1_est_NLR(:,index)-repmat(f1(index),K,1)).^2 ,2)))/sqrt(800)
std(sqrt(sum( (f1_est_NLR(:,index)-repmat(f1(index),K,1)).^2 ,2))/sqrt(800))

disp('mean and std f2 RMSE NLR:')
mean(sqrt(sum( (f2_est_NLR(:,index)-repmat(f2(index),K,1)).^2 ,2)))/sqrt(800)
std(sqrt(sum( (f2_est_NLR(:,index)-repmat(f2(index),K,1)).^2 ,2))/sqrt(800))

disp('mean and std f1 RMSE LR:')
mean(sqrt(sum( (f1_est_LR(:,index)-repmat(f1(index),K,1)).^2 ,2)))/sqrt(800)
std(sqrt(sum( (f1_est_LR(:,index)-repmat(f1(index),K,1)).^2 ,2))/sqrt(800))

disp('mean and std f2 RMSE LR:')
mean(sqrt(sum( (f2_est_LR(:,index)-repmat(f2(index),K,1)).^2 ,2)))/sqrt(800)
std(sqrt(sum( (f2_est_LR(:,index)-repmat(f2(index),K,1)).^2 ,2))/sqrt(800))

disp('mean and std f1 RMSE RDBR:')
mean(sqrt(sum( (f1_est_RDBR(:,index)-repmat(f1(index),K,1)).^2 ,2)))/sqrt(800)
std(sqrt(sum( (f1_est_RDBR(:,index)-repmat(f1(index),K,1)).^2 ,2))/sqrt(800))

disp('mean and std f2 RMSE RDBR:')
mean(sqrt(sum( (f2_est_RDBR(:,index)-repmat(f2(index),K,1)).^2 ,2)))/sqrt(800)
std(sqrt(sum( (f2_est_RDBR(:,index)-repmat(f2(index),K,1)).^2 ,2))/sqrt(800))

disp('mean and std f1 RMSE MMD:')
mean(sqrt(sum( (f1_est_MMD(:,index)-repmat(f1(index),K,1)).^2 ,2)))/sqrt(800)
std(sqrt(sum( (f1_est_MMD(:,index)-repmat(f1(index),K,1)).^2 ,2))/sqrt(800))

disp('mean and std f2 RMSE MMD:')
mean(sqrt(sum( (f2_est_MMD(:,index)-repmat(f2(index),K,1)).^2 ,2)))/sqrt(800)
std(sqrt(sum( (f2_est_MMD(:,index)-repmat(f2(index),K,1)).^2 ,2))/sqrt(800))

%---------------------------


vec1 = sqrt(sum(( f1_est_MMD(:,index)-repmat(f1(index),K,1)).^2 ,2))/sqrt(800);
vec2 = sqrt(sum(( f2_est_MMD(:,index)-repmat(f2(index),K,1)).^2 ,2))/sqrt(800);

f = f1+f2 + 0.3162*std(f1+f2)*randn(1,N);

% Graphics ----------------------------------------------------------
figure;subplot(6,2,1);
plot(t,f,'k')
text(0.005,5,'$f(t)$','interpreter','latex')

subplot(6,2,3)
plot(t,f1,'k');
ylim([-2 2.5])
text(0.005,1.8,'$s_1(t)$','interpreter','latex')
% 
subplot(6,2,4)
plot(t,f2,'k');
ylim([-1.5 6.5])
text(0.005,5.5,'$s_2(t)$','interpreter','latex')
% 
subplot(6,2,5)
fill([t fliplr(t)],[quantile(f1_est_NLR,0.025) fliplr(quantile(f1_est_NLR,0.975))],[0.8 0.8 0.8]); hold on;
plot(t,quantile(f1_est_NLR,0.025),'Color',[0.8 0.8 0.8]); plot(t,quantile(f1_est_NLR,0.975),'Color',[0.8 0.8 0.8])
plot(t,mean(f1_est_NLR),'k')
ylim([-2 2.5])
text(0.005,1.8,'$\tilde{s_1}(t)$ (SAMD)','interpreter','latex')

subplot(6,2,6)
fill([t fliplr(t)],[quantile(f2_est_NLR,0.025) fliplr(quantile(f2_est_NLR,0.975))],[0.8 0.8 0.8]); hold on;
plot(t,quantile(f2_est_NLR,0.025),'Color',[0.8 0.8 0.8]); plot(t,quantile(f2_est_NLR,0.975),'Color',[0.8 0.8 0.8])
plot(t,mean(f2_est_NLR),'k')
ylim([-1.5 6.5])
text(0.005,5.5,'$\tilde{s_2}(t)$ (SAMD)','interpreter','latex')

subplot(6,2,7)
fill([t fliplr(t)],[quantile(f1_est_LR,0.025) fliplr(quantile(f1_est_LR,0.975))],[0.8 0.8 0.8]); hold on;
plot(t,quantile(f1_est_LR,0.025),'Color',[0.8 0.8 0.8]); plot(t,quantile(f1_est_LR,0.975),'Color',[0.8 0.8 0.8])
plot(t,mean(f1_est_LR),'k')
ylim([-2 2.5])
text(0.005,1.8,'$\tilde{s_1}(t)$ (LR)','interpreter','latex')

subplot(6,2,8)
fill([t fliplr(t)],[quantile(f2_est_LR,0.025) fliplr(quantile(f2_est_LR,0.975))],[0.8 0.8 0.8]); hold on;
plot(t,quantile(f2_est_LR,0.025),'Color',[0.8 0.8 0.8]); plot(t,quantile(f2_est_LR,0.975),'Color',[0.8 0.8 0.8])
plot(t,mean(f2_est_LR),'k')
ylim([-1.5 6.5])
text(0.005,5.5,'$\tilde{s_2}(t)$ (LR)','interpreter','latex')

subplot(6,2,9)
fill([t fliplr(t)],[quantile(f1_est_RDBR,0.025) fliplr(quantile(f1_est_RDBR,0.975))],[0.8 0.8 0.8]); hold on;
plot(t,quantile(f1_est_RDBR,0.025),'Color',[0.8 0.8 0.8]); plot(t,quantile(f1_est_RDBR,0.975),'Color',[0.8 0.8 0.8])
plot(t,mean(f1_est_RDBR),'k')
ylim([-2 2.5])
text(0.005,1.8,'$\tilde{s_1}(t)$ (RDBR)','interpreter','latex')

subplot(6,2,10)
fill([t fliplr(t)],[quantile(f2_est_RDBR,0.025) fliplr(quantile(f2_est_RDBR,0.975))],[0.8 0.8 0.8]); hold on;
plot(t,quantile(f2_est_RDBR,0.025),'Color',[0.8 0.8 0.8]); plot(t,quantile(f2_est_RDBR,0.975),'Color',[0.8 0.8 0.8])
plot(t,mean(f2_est_RDBR),'k')
ylim([-1.5 6.5])
text(0.005,5.5,'$\tilde{s_2}(t)$ (RDBR)','interpreter','latex')

subplot(6,2,11)
fill([t fliplr(t)],[quantile(f1_est_MMD(vec2<10,:),0.025) fliplr(quantile(f1_est_MMD(vec2<10,:),0.975))],[0.8 0.8 0.8]); hold on;
plot(t,quantile(f1_est_MMD(vec2<10,:),0.025),'Color',[0.8 0.8 0.8]); plot(t,quantile(f1_est_MMD(vec2<10,:),0.975),'Color',[0.8 0.8 0.8])
plot(t,mean(f1_est_MMD(vec2<10,:)),'k')
ylim([-2 2.5])
text(0.005,1.8,'$\tilde{s_1}(t)$ (MMD)','interpreter','latex')
xlabel('$t$','interpreter','latex')

subplot(6,2,12)
fill([t fliplr(t)],[quantile(f2_est_MMD(vec2<10,:),0.025) fliplr(quantile(f2_est_MMD(vec2<10,:),0.975))],[0.8 0.8 0.8]); hold on;
plot(t,quantile(f2_est_MMD(vec2<10,:),0.025),'Color',[0.8 0.8 0.8]); plot(t,quantile(f2_est_MMD(vec2<10,:),0.975),'Color',[0.8 0.8 0.8])
plot(t,mean(f2_est_MMD(vec2<10,:)),'k')
ylim([-1.5 6.5])
text(0.005,5.5,'$\tilde{s_2}(t)$ (MMD)','interpreter','latex')
xlabel('$t$','interpreter','latex')

