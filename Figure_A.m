% This code illustrates the usage of SAMD and LR on the "first simulated
% signal" from the paper "Signal decomposition for time-varying wave-shape 
% "functions and its biomedical application" by Marcelo A. Colominas and
% Hau-Tieng Wu.
%
% Marcelo A. Colominas
% email: macolominas@conicet.gov.ar
% 17-SEP-2020

t = 0:0.001:1-0.001; %temporal variable
N = length(t);

% first component
phi1 = 2*pi*20*t + 2*cos(2*2*pi*t);
ti = diff( mod  ( 0.5*phi1/pi,1 )  );
ti = find(ti<0);
f1 = zeros(size(t));
for i = 1:length(ti)
    f1 = f1 + exp(-200000*(t-t(ti(i))).^2);
end
f1 = f1 - mean(f1);

% second component
phi2 = (2*pi*10*t + 2*pi*5*t.^2);
f2 = 0.5*(cos(phi2) + 0.75*cos(2.05*phi2));

f = f1 + f2;% you can try adding noise here

D = [10 2];

order_Phi = 3;

tic;
[f_est,modos_est,b_e] = SAMD(f,[ones(size(t));ones(size(t))],[phi1;phi2],D,order_Phi);
toc

tic;
[f_est_lin,modos_est_lin] = LR(f,[ones(size(t));ones(size(t))],[phi1;phi2],D);
toc

%---Errors------------------------

norm(f1-modos_est(1,:))/sqrt(N)
norm(f2-modos_est(2,:))/sqrt(N)

norm(f1-modos_est_lin(1,:))/sqrt(N)
norm(f2-modos_est_lin(2,:))/sqrt(N)



%---Graphics----------------------

figure;subplot(6,2,1);
plot(t,f,'k')
text(0.005,1.5,'$f(t)$','interpreter','latex')

subplot(6,2,3)
plot(t,f1,'k');
ylim([-0.08 1.25])
text(0.005,1.05,'$s_1(t)$','interpreter','latex')

subplot(6,2,4)
plot(t,f2,'k');
ylim([-0.9 1.35])
text(0.005,1.05,'$s_2(t)$','interpreter','latex')

subplot(6,2,5)
plot(t,f1,'r-.');hold on; plot(t,modos_est(1,:),'k')
ylim([-0.1 1.25])
text(0.005,1.05,'$\tilde{s_1}(t)$ (SAMD)','interpreter','latex')

subplot(6,2,6)
plot(t,f2,'r-.');hold on; plot(t,modos_est(2,:),'k')
ylim([-0.9 1.35])
text(0.005,1.05,'$\tilde{s_2}(t)$ (SAMD)','interpreter','latex')

subplot(6,2,7)
plot(t,f1,'r-.');hold on; plot(t,modos_est_lin(1,:),'k')
ylim([-0.17 1.25])
text(0.005,1.05,'$\tilde{s_1}(t)$ (LR)','interpreter','latex')

subplot(6,2,8)
plot(t,f2,'r-.');hold on; plot(t,modos_est_lin(2,:),'k')
ylim([-0.9 1.35])
text(0.005,1.05,'$\tilde{s_2}(t)$ (LR)','interpreter','latex')

