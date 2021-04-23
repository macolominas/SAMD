% This code produces Fig. 9 of the Supplemental Material of the paper "Decomposing non-stationary signals
% with time-varying wave-shape functions" by Marcelo A. Colominas and Hau-Tieng Wu.
%
%
% Marcelo A. Colominas
% email: macolominas@conicet.gov.ar
% 23-APR-2021
addpath(genpath('DeCom-master')); %to compare with RDBR and MMD (developed by Haizhao Yang)

t = 0:0.001:1-0.001;%temporal variable
N = length(t);

phi1 = 2*(2*pi*3*t + 2*pi*3*t.^2);%phase function
phi2 = 2*(2*pi*5*t + 2*pi*3.5*t.^2 + 0.5*cos(2*pi*t));%phase function

% Constructing s_1
f1 = 1.5*cos(phi1) + 0.25*cos(2.05*phi1) + 0.1*cos(3.05*phi1 + 0.01*phi1.^2) + 0.1*cos(4.05*phi1 + 0.01*phi1.^2); 
for i = 5:10
    f1 = f1 + 0.1*cos((i+0.05)*phi1+ 0.01*phi1.^2);
end;

% Constructing s_2
f2 = cos(phi2);
for i = 2:10
    f2 = f2 + sqrt((1/i))*cos((i+0.01)*phi2);
end;

f = f1 + f2;

%-------Estimating ridges and phases---------------------------------

bt = 1:N;
redun = 4;
gamma = 0;
sigma = 0.25;
fmax = 0.1;
ft = 1/redun:1/redun:round(fmax*N);

[STFT,SST1,SST2,SST3,SST4,omega,omega2,omega3,omega4,tau2,tau3,phi22p,phi23p,phi24p] = sstn_test_modL_new(f,gamma,sigma,ft,bt,redun);
alpha = 3;
(1/(1-alpha))*log(sum(sum((abs(STFT./sum(abs(STFT(:))))).^alpha)))

RTF = SST2;
jump = redun/2;
d = 4;
c1 = exridge(RTF,0,0,jump);
plot(c1,'r')

for i = 1:N
    aux(i) = sum(RTF(c1(i)-d:c1(i)+d,i));
    RTF(c1(i)-d:c1(i)+d,i) = 0;
end;
RTF(1:10,:) = 0;
phi1_est = phase(aux);
A1_est = ones(size(t));

c2 = exridge(RTF,0,0,jump);
plot(c1); plot(c2)
for i = 1:N
    aux(i) = sum(RTF(c2(i)-d:c2(i)+d,i));
end;
phi2_est = phase(aux);
A2_est = ones(size(t));
%--------------------------------------------------------------------

%--Parameters for our methods----------------------------------------
D1 = 10;
D2 = 10;
order_Phi = 3;
%--------------------------------------------------------------------

tic;
warning('off')
[f_est,f1_est,f2_est,b_e] = two_wsf_nonlinear(f,A1_est,phi1_est,D1,A2_est,phi2_est,D2,order_Phi);
warning('on')
toc

tic;
[f_est_lin,f1_est_lin,f2_est_lin,~,~,coef_est] = two_wsf_linear(f,A1_est,phi1_est,D1,A2_est,phi2_est,D2);
toc

%---RDBR-------------------------------------------------------------

N = length(t);

opt.maxiter = 3000;
opt.eps_error = 1e-12;
opt.show = 0;
opt.nknots = 20;
opt.knotremoval_factor= 1.01;
opt.order = 3;
opt.eps_diff = opt.eps_error;
    
ins_pre_phase = [phi1_est;phi2_est]*0.5/pi;
ins_amplt = [A1_est;A2_est];
    
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
    [shapeInLoop,comp,errorRec,SL2Rec,iter,flag] = srcIterRegJC(f,N,numGroup,ins_amplt,ins_pre_phase,opt,fTrue,shapeTrue);
    toc
end
        
%--------------------------------------------------------------------

%---MMD--------------------------------------------------------------
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
[shape,component,Hcoef,flag,idx,iter,iterDR] = DeCom_MMD(f,t,2,[A1_est;A2_est],inst_freq,[phi1_est;phi2_est]*0.5/pi,opt);
toc

%--------------------------------------------------------------------
index = 101:900;

norm(f1(index)-f1_est(index))/sqrt(800)
norm(f2(index)-f2_est(index))/sqrt(800)
norm(f1(index)-f1_est_lin(index))/sqrt(800)
norm(f2(index)-f2_est_lin(index))/sqrt(800)
norm(f1(index)-comp{1}(index))/sqrt(800)
norm(f2(index)-comp{2}(index))/sqrt(800)
norm(f1(index)-component{1}(index))/sqrt(800)
norm(f2(index)-component{2}(index))/sqrt(800)

% ----- Graphics ----------------------------------------------------
figure;subplot(6,2,1);
plot(t,f,'k')
text(0.005,5,'$f(t)$','interpreter','latex')

subplot(6,2,3)
plot(t,f1,'k');
ylim([-2 2.5])
text(0.005,1.8,'$s_1(t)$','interpreter','latex')

subplot(6,2,4)
plot(t,f2,'k');
ylim([-1.5 6.5])
text(0.005,5.5,'$s_2(t)$','interpreter','latex')

subplot(6,2,5)
plot(t,f1,'r-.');hold on; plot(t,f1_est,'k')
ylim([-2 2.5])
text(0.005,1.8,'$\tilde{s_1}(t)$ (SAMD)','interpreter','latex')

subplot(6,2,6)
plot(t,f2,'r-.');hold on; plot(t,f2_est,'k')
ylim([-1.5 6.5])
text(0.005,5.5,'$\tilde{s_2}(t)$ (SAMD)','interpreter','latex')

subplot(6,2,7)
plot(t,f1,'r-.');hold on; plot(t,f1_est_lin,'k')
ylim([-2 2.5])
text(0.005,1.8,'$\tilde{s_1}(t)$ (LR)','interpreter','latex')

subplot(6,2,8)
plot(t,f2,'r-.');hold on; plot(t,f2_est_lin,'k')
ylim([-1.5 6.5])
text(0.005,5.5,'$\tilde{s_2}(t)$ (LR)','interpreter','latex')

subplot(6,2,9)
plot(t,f1,'r-.');hold on; plot(t,comp{1},'k')
ylim([-2 2.5])
text(0.005,1.8,'$\tilde{s_1}(t)$ (RDBR)','interpreter','latex')

subplot(6,2,10)
plot(t,f2,'r-.');hold on; plot(t,comp{2},'k')
ylim([-1.5 6.5])
text(0.005,5.5,'$\tilde{s_2}(t)$ (RDBR)','interpreter','latex')

subplot(6,2,11)
plot(t,f1,'r-.');hold on; plot(t,component{1},'k')
xlabel('$t$','interpreter','latex')
ylim([-2 2.5])
text(0.005,1.8,'$\tilde{s_1}(t)$ (MMD)','interpreter','latex')

subplot(6,2,12)
plot(t,f2,'r-.');hold on; plot(t,component{2},'k')
xlabel('$t$','interpreter','latex')
ylim([-1.5 6.5])
text(0.005,5.5,'$\tilde{s_2}(t)$ (MMD)','interpreter','latex')


