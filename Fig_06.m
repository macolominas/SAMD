% This code produces Fig. 6 from the paper "Decomposing non-stationary signals
% with time-varying wave-shape functions" by Marcelo A. Colominas and Hau-Tieng Wu.
%
% The running of this file might take several hours. Please, modify
% variable I (line 70) accordingly.
%
% Marcelo A. Colominas
% email: macolominas@conicet.gov.ar
% 23-APR-2021
addpath(genpath('DeCom-master')); %to compare with RDBR and MMD (developed by Haizhao Yang)

t = 0:0.001:1-0.001;%temporal variable
N = length(t);

phi1 = 2*(2*pi*3*t + 2*pi*3*t.^2); %phase function
phi2 = 2*(2*pi*5*t + 2*pi*3.5*t.^2 + 0.5*cos(2*pi*t)); %phase function

% Parameters for ridges, amplitudes and phases estimations-----------
bt = 1:N;
redun = 4;
gamma = 0;
sigma = 0.25;
fmax = 0.05;
ft = 1/redun:1/redun:round(fmax*N);

index = 101:900;
jump = redun;
d = 4;
%--------------------------------------------------------------------

%--Parameters for our methods----------------------------------------
D1 = 10;
D2 = 10;
order_Phi = 3;
%--------------------------------------------------------------------

%--Parameters for MMD------------------------------------------------
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
%--------------------------------------------------------------------

%--Parameters for RDBR ----------------------------------------------
optR.maxiter = 3000;
optR.eps_error = 1e-12;
optR.show = 0;
optR.nknots = 20;
optR.knotremoval_factor= 1.01;
optR.order = 3;
optR.eps_diff = opt.eps_error;

%--------------------------------------------------------------------

I = 100;%number of realizations

error_1_LR = zeros(3,I);
error_1_SAMD = zeros(3,I);
error_1_RDBR = zeros(3,I);
error_1_MMD = zeros(3,I);

error_2_LR = zeros(3,I);
error_2_SAMD = zeros(3,I);
error_2_RDBR = zeros(3,I);
error_2_MMD = zeros(3,I);

tiempo_LR = zeros(3,I);
tiempo_SAMD = zeros(3,I);
tiempo_RDBR = zeros(3,I);
tiempo_MMD = zeros(3,I);

parfor i=1:I

noise = zeros(3,N);
    
shift_f1 = 0.05*rand(1,9)+0.05;% variables \xi_p, for p = 2,...,10.

% Construction of Y(t)
ff = cumsum(randn(1,N));
IFr = ff./max(abs(ff));
IF = smooth(IFr,80);

phi1_ = (4*pi*3*t + 4*pi*3*t.^2) + 2*pi*(cumsum(IF)')/1000;%phase function
%-----
% Construction of s_1
f1 = 1.5*cos(phi1_) + 0.25*cos((2+shift_f1(1))*phi1_) + 0.1*cos((3+shift_f1(2))*phi1_ + 0.01*phi1_.^2) + 0.1*cos((4+shift_f1(3))*phi1_ + 0.01*phi1_.^2);
for p = 5:10
    f1 = f1 + 0.1*cos((p+0.05)*phi1_+ 0.01*phi1_.^2);
end;

% Construction of Z(t)
ff = cumsum(randn(1,N));
IFr = ff./max(abs(ff));
IF = smooth(IFr,50);
phi2_ = phi1_ + 8*pi*t + 2*pi*t.^2 + 0.5*cos(2*pi*t) + 2*pi*(cumsum(IF)')/1000;%phase function

% Construction of s_2
f2 = cos(phi2_);
shift_f2 = 0.01*rand(1,9)+0.01;% variables \zeta_p, for p = 2,...,10.
for p = 2:10
    f2 = f2 + sqrt((1/p))*cos((p+shift_f2(p-1))*phi2_);
end;
%--------------------------------------------------------------------

noise(1,:) = randn(1,N);% different noises

% ARMA Noise --------------------------------------------------------
e = random('t',4,size(t)) ; e = e(:) ; % Student T-4
aux = zeros(1,N);
aux(1) = e(1);
for k = 2:N
    aux(k) = 0.5*aux(k-1) + e(k) + 0.5*e(k-1);
end;
aux = aux/std(aux);
noise(2,:) = aux;%ARMA
% -------------------------------------------------------------------

aux = poissrnd(1,[1,N]);
aux = aux/std(aux);
noise(3,:) = aux;%Poisson noise

for j = 1:3

f = f1 + f2 + 0.3162*std(f1+f2)*noise(j,:);

[STFT,SST1,SST2,SST3,SST4,omega,omega2,omega3,omega4,tau2,tau3,phi22p,phi23p,phi24p] = sstn_test_modL_new(f,gamma,sigma,ft,bt,redun);

RTF = SST2;
c1 = exridge(RTF,0,0,jump);

for n = 1:N
    a = max(1,c1(n)-d);
    b = min(200,c1(n)+d);
    aux(n) = sum(RTF(a:b,n));
    RTF(a:b,n) = 0;
end;
RTF(1:10,:) = 0;
phi1_est = phase(aux);
A1_est = abs(aux);
c2 = exridge(RTF,0,0,jump);

for n = 1:N
    a = max(1,c2(n)-d);
    b = min(200,c2(n)+d);
    aux(n) = sum(RTF(a:b,n));
end;
phi2_est = phase(aux);
A2_est = abs(aux);

%--------------------------------------------------------------------


tic;
warning('off')
[f_est,f1_est,f2_est,b_e] = two_wsf_nonlinear(f,ones(size(t)),phi1_est,D1,ones(size(t)),phi2_est,D2,order_Phi);
warning('on')%OJO ACAAAAAAAAAAAAAAAAAAAAAAAAA^AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
tiempo_SAMD(j,i) = toc

tic;
[f_est_lin,f1_est_lin,f2_est_lin,~,~,coef_est] = two_wsf_linear(f,ones(size(t)),phi1_est,D1,ones(size(t)),phi2_est,D2);
tiempo_LR(j,i) = toc

%---RDBR-------------------------------------------------------------
    
ins_pre_phase = [phi1_est;phi2_est]*0.5/pi;
ins_amplt = [ones(size(t));ones(size(t))];
    
if (1) 
    fTrue = cell(1,2);
    fTrue{1} = f1;
    fTrue{2} = f2;
    numGroup = 2;
    fff = f1 + f2;
    shapeTrue = cell(1,2);
%         shapeTrue{1} = @(x) sh1(x);
%         shapeTrue{2} = @(x) sh2(x);
    tic;
    [shapeInLoop,comp,errorRec,SL2Rec,iter,flag] = srcIterRegJC(f,N,numGroup,ins_amplt,ins_pre_phase,optR,fTrue,shapeTrue);
    tiempo_RDBR(j,i) = toc
end    
%--------------------------------------------------------------------

% %---MMD------------------------------------------------------------

    
tic;
inst_freq = zeros(2,N);
inst_freq(1,1:end-1) = N*diff(phi1_est)*0.5/pi; inst_freq(1,end) = inst_freq(1,end-1);
inst_freq(2,1:end-1) = N*diff(phi2_est)*0.5/pi; inst_freq(2,end) = inst_freq(2,end-1);
[shape,component,Hcoef,flag,idx,iter,iterDR] = DeCom_MMD(f,t,2,[ones(size(t));ones(size(t))],inst_freq,[phi1_est;phi2_est]*0.5/pi,opt);
tiempo_MMD(j,i) = toc

%--------------------------------------------------------------------

error_1_SAMD(j,i) = norm(f1(index)-f1_est(index))/sqrt(800);
error_2_SAMD(j,i) = norm(f2(index)-f2_est(index))/sqrt(800);
error_1_LR(j,i) = norm(f1(index)-f1_est_lin(index))/sqrt(800);
error_2_LR(j,i) = norm(f2(index)-f2_est_lin(index))/sqrt(800);
error_1_RDBR(j,i) = norm(f1(index)-comp{1}(index))/sqrt(800);
error_2_RDBR(j,i) = norm(f2(index)-comp{2}(index))/sqrt(800);
error_1_MMD(j,i) = norm(f1(index)-component{1}(index))/sqrt(800);
error_2_MMD(j,i) = norm(f2(index)-component{2}(index))/sqrt(800);

end;%j
i
end;%i

% Graphics
figure;
subplot(331)
boxplot([error_1_SAMD(1,:)' error_1_LR(1,:)' error_1_RDBR(1,:)' error_1_MMD(1,:)'],'Labels',{'SAMD','LR','RDBR','MMD'})
hold on
ylim([0.2 1.5])
plot([1 2],[1.2 1.2],'k')
text(1.4,1.215,'*')
plot([1 3],[1.3 1.3],'k')
text(1.9,1.315,'*')
plot([1 4],[1.4 1.4],'k')
text(2.4,1.415,'*')
title('Gaussian noise (10 dB)')
ylabel('RMSE for $s_1(t)$','interpreter','latex')

subplot(332)
boxplot([error_1_SAMD(2,:)' error_1_LR(2,:)' error_1_RDBR(2,:)' error_1_MMD(2,:)'],'Labels',{'SAMD','LR','RDBR','MMD'})
hold on
ylim([0.2 1.5])
plot([1 4],[1.4 1.4],'k')
text(2.4,1.415,'*')
title('ARMA(1,1) noise (10 dB)')

subplot(333)
boxplot([error_1_SAMD(3,:)' error_1_LR(3,:)' error_1_RDBR(3,:)' error_1_MMD(3,:)'],'Labels',{'SAMD','LR','RDBR','MMD'})
hold on
ylim([0.2 1.5])
plot([1 4],[1.4 1.4],'k')
text(2.4,1.415,'*')
title('Poisson noise (10 dB)')

subplot(334)
boxplot([error_2_SAMD(1,:)' error_2_LR(1,:)' error_2_RDBR(1,:)' error_2_MMD(1,:)'],'Labels',{'SAMD','LR','RDBR','MMD'})
hold on
plot([1 2],[1.2 1.2],'k')
text(1.4,1.215,'*')
plot([1 3],[1.3 1.3],'k')
text(1.9,1.315,'*')
plot([1 4],[1.4 1.4],'k')
text(2.4,1.415,'*')
ylim([0.2 1.5])
ylabel('RMSE for $s_2(t)$','interpreter','latex')

subplot(335)
boxplot([error_2_SAMD(2,:)' error_2_LR(2,:)' error_2_RDBR(2,:)' error_2_MMD(2,:)'],'Labels',{'SAMD','LR','RDBR','MMD'})
hold on
plot([1 2],[1.2 1.2],'k')
text(1.4,1.215,'*')
plot([1 3],[1.3 1.3],'k')
text(1.9,1.315,'*')
ylim([0.2 1.5])

subplot(336)
boxplot([error_2_SAMD(3,:)' error_2_LR(3,:)' error_2_RDBR(3,:)' error_2_MMD(3,:)'],'Labels',{'SAMD','LR','RDBR','MMD'})
hold on
plot([1 2],[1.2 1.2],'k')
text(1.4,1.215,'*')
plot([1 3],[1.3 1.3],'k')
text(1.9,1.315,'*')
plot([1 4],[1.4 1.4],'k')
text(2.4,1.415,'*')
ylim([0.2 1.5])


