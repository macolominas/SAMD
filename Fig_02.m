% This code produces Fig. 3 from the paper "Decomposing non-stationary signals
% with time-varying wave-shape functions" 
% by Marcelo A. Colominas and Hau-Tieng Wu.
%
%
% Marcelo A. Colominas
% email: macolominas@conicet.gov.ar
% 23-APR-2021
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

%-------Parameters for Estimating ridges and phases-----------------
bt = 1:N;
redun = 4;
gamma = 0;
sigma = 0.25;
fmax = 0.1;
ft = 1/redun:1/redun:round(fmax*N);
%--------------------------------------------------------------------

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


% Graphics
figure;
subplot(331)
imagesc(t,ft,(abs(STFT))); colormap(1-gray)
ylabel('frequency [Hz]')
xlabel('time')
title('|STFT|')

subplot(332)
imagesc(t,ft,(abs(SST2))); colormap(1-gray); hold on
xlabel('time')
title('|SST2|')

subplot(333)
imagesc(t,ft(1:30*redun),(abs(SST2(1:30*redun,:)))); colormap(1-gray); hold on
plot(t(1:end-1),N*diff(phi1)/(2*pi),'b')
plot(t,c1/redun,'b--','Linewidth',2)
plot(t(1:end-1),N*diff(phi2)/(2*pi),'r')
plot(t,c2/redun,'r--','Linewidth',2)
xlabel('time')
title('|SST2| with ridges')