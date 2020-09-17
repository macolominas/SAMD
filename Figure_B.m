% This code illustrates the usage of SAMD and LR on the "second simulated
% signal" from the paper "Signal decomposition for time-varying wave-shape 
% "functions and its biomedical application" by Marcelo A. Colominas and
% Hau-Tieng Wu.
%
% Marcelo A. Colominas
% email: macolominas@conicet.gov.ar
% 17-SEP-2020

t = 0:0.001:1-0.001; %temporal variable

N = length(t); 

% third component
phi3 = (2*pi*18*t + 4*pi*t.^2 + cos(2*2*pi*t)); 
ti = diff( mod  ( 0.5*phi3/pi,1 )  );
ti = find(ti<0);
M = 200000;
f3 = zeros(size(t));
for i = 1:length(ti)
    f3 = f3 - 0.5*2*2000*(t-t(ti(i))).*exp(-M*(t-t(ti(i))).^2);
end;
f3 = f3 - mean(f3);

% second component
phi2 = 2*pi*(14*t+t.^2+t.^3);
f2 = cos(phi2) + 0.75*cos(2.05*phi2) + 0.25*cos(2.95*phi2);

% first component
phi1 = 2*pi*(12*t+2*t.^2);
f1 = cos(phi1) + 0.5*cos(1.95*phi1 + 0.0001*phi1.^2);

% mixed signal
f = f1 + f2 + f3;% you can try adding noise here

D = [2 3 20];

order_Phi = 3;

tic;
[f_est,modos_est,b_e] = SAMD(f,[ones(size(t));ones(size(t));ones(size(t))],[phi1;phi2;phi3],D,order_Phi);
toc

tic;
[f_est_lin,modos_est_lin] = LR(f,[ones(size(t));ones(size(t));ones(size(t))],[phi1;phi2;phi3],D);
toc

%---Errors------------------------

norm(f1-modos_est(1,:))/sqrt(N)
norm(f2-modos_est(2,:))/sqrt(N)
norm(f3-modos_est(3,:))/sqrt(N)

norm(f1-modos_est_lin(1,:))/sqrt(N)
norm(f2-modos_est_lin(2,:))/sqrt(N)
norm(f3-modos_est_lin(3,:))/sqrt(N)



%---Graphics----------------------
figure;subplot(6,3,1);
plot(t,f,'k')
text(0.005,4,'$f(t)$','interpreter','latex')

subplot(6,3,4)
plot(t,f1,'k');
% ylim([-0.08 1.25])
text(0.005,1.5,'$s_1(t)$','interpreter','latex')

subplot(6,3,5)
plot(t,f2,'k');
ylim([-2 3])
text(0.005,2.5,'$s_2(t)$','interpreter','latex')

subplot(6,3,6)
plot(t,f3,'k');
ylim([-2 3])
text(0.005,2.5,'$s_3(t)$','interpreter','latex')

subplot(6,3,7)
plot(t,f1,'r-.');hold on; plot(t,modos_est(1,:),'k')
% ylim([-0.1 1.25])
text(0.005,1.5,'$\tilde{s_1}(t)$ (SAMD)','interpreter','latex')

subplot(6,3,8)
plot(t,f2,'r-.');hold on; plot(t,modos_est(2,:),'k')
ylim([-2 3])
text(0.005,2.5,'$\tilde{s_2}(t)$ (SAMD)','interpreter','latex')

subplot(6,3,9)
plot(t,f3,'r-.');hold on; plot(t,modos_est(3,:),'k')
ylim([-2 3])
text(0.005,2.5,'$\tilde{s_2}(t)$ (SAMD)','interpreter','latex')

subplot(6,3,10)
plot(t,f1,'r-.');hold on; plot(t,modos_est_lin(1,:),'k')
%ylim([-0.17 1.25])
text(0.005,1.5,'$\tilde{s_1}(t)$ (LR)','interpreter','latex')
 
subplot(6,3,11)
plot(t,f2,'r-.');hold on; plot(t,modos_est_lin(2,:),'k')
ylim([-2 3])
text(0.005,2.5,'$\tilde{s_2}(t)$ (LR)','interpreter','latex')

subplot(6,3,12)
plot(t,f3,'r-.');hold on; plot(t,modos_est_lin(3,:),'k')
ylim([-2 3])
text(0.005,2.5,'$\tilde{s_2}(t)$ (LR)','interpreter','latex')

