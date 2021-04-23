function [STFT,SST1,SST2,SST3,SST4,omega,omega2,omega3,omega4,tau2,tau3,phi22p,phi23p,phi24p] = sstn_test_modL(s,gamma,sigma,ft,bt,redun)
% sstn : computes the STFT of a signal and different versions of synchrosqueezing/reassignment.
%   Uses a Gaussian window.
%
% INPUTS:
%   s: real or complex signal, must be of length 2^N
%   gamma: threshold
%   sigma: window parameter
%   ft: frequency bins
%   bt: time bins
%
% OUTPUTS:
%   STFT: the short-time Fourier transform
%   SST1: standard synchrosqueezing
%   SST2: vertical second-order synchrosqueezing
%   SST3: vertical third-order synchrosqueezing
%   SST4: vertical fourth-order synchrosqueezing
%   omega: instantaneous frequency (vertical reassignment operator)
%   tau2: second-order phase delay (horizontal reassignment operator)
%   tau3: third-order phase delay (horizontal reassignment operator)
%   omega2: second-order instantaneous frequency
%   omega3: third-order instantaneous frequency

% Consume menos memoria porque calcula fila a fila todas las matrices
% involucradas a medida que son necesarias. Al no guardar todas las STFTs,
% y calcular sólo la fila a utilizar en esa iteración, utiliza menos memoria
% y permite analizar señales más largas (2^17 muestras por ejemplo).

% checking length of signal
% n = length(s);
% nv = log2(n);
% if mod(nv,1)~=0
%     warning('The signal is not a power of two, truncation to the next power');
%     s = s(1:2^floor(nv));
% end
n = length(s);
%s = s(:);

% Optional parameters
if nargin<5
    ft = 1:n/2;
    bt = 1:n;
end
nb = length(bt);
neta = length(ft);
%sz=zeros(n,1);
% Padding
% sleft = flipud(conj(sz(2:n/2+1)));
% sright = flipud(sz(end-n/2:end-1));
% x = [sleft; s ; sright];
% clear xleft xright;

% Window definition
% t = -0.5:1/n:0.5-1/n;t=t';
% g =  1/sigma*exp(-pi/sigma^2*t.^2);
% gp = -2*pi/sigma^2*t .* g; % g'
%gpp = (-2*pi/sigma^2+4*pi^2/sigma^4*t.^2) .* g; % g''
S = fftshift(fft(s));
f = -n/2:n/2-1;

% Initialization
STFT = zeros(neta,nb);
SST1 = zeros(neta,nb);
SST2 = zeros(neta,nb);
SST3 = zeros(neta,nb);
SST4 = zeros(neta,nb);
omega = zeros(neta,nb);
tau2 = zeros(neta,nb);
tau3 = zeros(neta,nb);
tau4 = zeros(neta,nb);
omega2 = zeros(neta,nb);
omega3 = zeros(neta,nb);
omega4 = zeros(neta,nb);
phi22p = complex(zeros(neta,nb));
%phi2p = zeros(neta,nb);
phi23p = complex(zeros(neta,nb));
phi33p = complex(zeros(neta,nb));
%phi3p = zeros(neta,nb);
phi24p = complex(zeros(neta,nb));
phi34p = complex(zeros(neta,nb));
phi44p = complex(zeros(neta,nb));

vg = complex(zeros(8,nb)); % Este consume menos memoria.
vgp = complex(zeros(6,nb));
Y = complex(zeros(7,nb,7));



for b = 1:length(ft) % For each row
    
    %% Computes STFT and reassignment operators
    for i = 0:7
        g_hat = (sigma*sqrt(pi)/(1i*2*pi))^i*...
            hermite_poly(sigma*sqrt(pi)*(f-ft(b)),i).*...
            exp(-sigma^2*pi*(f-ft(b)).^2);
        %          plot(f,abs(S)); hold on; plot(f,abs(g_hat));
        %          clf;
         tmp = ifft(S.*conj(g_hat));
         vg(i+1,:) = tmp(bt);
    end
    
    for i = 0:5
        tgp_hat = (-2*pi/sigma^2)*(sigma*sqrt(pi)/(1i*2*pi))^(i+1)*...
            hermite_poly(sigma*sqrt(pi)*(f-ft(b)),i+1).*...
            exp(-sigma^2*pi*(f-ft(b)).^2);
        tmp = ifft(S.*conj(tgp_hat));
        vgp(i+1,:) = tmp(bt);
    end
    
    
    %% second-order operator tau
    tau2(b,:) = vg(2,:)./vg(1,:);
    % third order operator tau
    tau3(b,:) = vg(3,:)./vg(1,:);
    % four order operator tau
    tau4(b,:) = vg(4,:)./vg(1,:);
    
    %% Y expressions
    for i = 1:7
        for j = 1:7
            if i>=j
                Y(i,:,j) = vg(1,:).*vg(i+1,:) - vg(j,:).*vg(i-j+2,:);
            end
        end
    end
    

    %% W expressions
    W2 = 1/2/1i/pi*(vg(1,:).^2+vg(1,:).*vgp(2,:)-vg(2,:).*vgp(1,:));
    W3 = 1/2/1i/pi*(2*vg(1,:).*vg(2,:)+vg(1,:).*vgp(3,:)-vg(3,:).*vgp(1,:));
    W4 = 1/2/1i/pi*(2*vg(1,:).*vg(3,:)+2*vg(2,:).^2+vg(1,:).*vgp(4,:) - vg(4,:).*vgp(1,:)+vg(2,:).*vgp(3,:) - vg(3,:).*vgp(2,:));
    
    %% operator omega
    omega(b,:) = (ft(b))-real(vgp(1,:)/2/1i/pi./vg(1,:));%ac� rest� 1 (MAC)
    
    %% operator hat p: estimations of frequency modulation
    %SST2
    phi22p(b,:) = W2./Y(2,:,2);
    omega2(b,:) = omega(b,:) + real(phi22p(b,:).*tau2(b,:));
    
    %SST3
    phi33p(b,:) = (W3.*Y(2,:,2)-W2.*Y(3,:,3))./(Y(4,:,3).*Y(2,:,2)-Y(3,:,2).*Y(3,:,3));
    phi23p(b,:) = W2./Y(2,:,2) - phi33p(b,:).*Y(3,:,2)./Y(2,:,2);
    omega3(b,:) = omega(b,:) + real(phi23p(b,:).*tau2(b,:))+ real(phi33p(b,:).*tau3(b,:));
    
    %SST4
    phi44p(b,:) =((Y(4,:,3).*Y(2,:,2)-Y(3,:,2).*Y(3,:,3)).*W4-(W3.*Y(2,:,2)-W2.*Y(3,:,3)).*(Y(5,:,4)+Y(5,:,3)-Y(5,:,2))+(W3.*Y(3,:,2)-W2.*Y(4,:,3)).*(Y(4,:,4)+Y(4,:,3)-Y(4,:,2)))...
        ./((Y(4,:,3).*Y(2,:,2)-Y(3,:,2).*Y(3,:,3)).*(Y(6,:,4)+Y(6,:,3)-Y(6,:,2))-(Y(5,:,3).*Y(2,:,2)-Y(4,:,2).*Y(3,:,3)).*(Y(5,:,4)+Y(5,:,3)-Y(5,:,2))+(Y(5,:,3).*Y(3,:,2)-Y(4,:,2).*Y(4,:,3)).*(Y(4,:,4)+Y(4,:,3)-Y(4,:,2)));
    
    phi34p(b,:) = (W3.*Y(2,:,2)-W2.*Y(3,:,3))./(Y(4,:,3).*Y(2,:,2)-Y(3,:,2).*Y(3,:,3))-phi44p(b,:).*(Y(5,:,3).*Y(2,:,2)-Y(4,:,2).*Y(3,:,3))./(Y(4,:,3).*Y(2,:,2)-Y(3,:,2).*Y(3,:,3));
    phi24p(b,:) = W2./Y(2,:,2) - phi34p(b,:).*Y(3,:,2)./Y(2,:,2) - phi44p(b,:).*Y(4,:,2)./Y(2,:,2);
    omega4(b,:) = omega(b,:) + real(phi24p(b,:).*tau2(b,:))+ real(phi34p(b,:).*tau3(b,:))+ real(phi44p(b,:).*tau4(b,:));
    
    % Storing STFT
    STFT(b,:) = vg(1,:).*exp(1i*pi*(bt-1)); % compensates the tranlation 1/2 of s
    
end

%% Reassignment step
for b=1:nb
    for eta=1:neta
        if abs(STFT(eta,b))>gamma
           
%%%%%%%%%%%%%%%%%%%%%%%%%%SST1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            k = 1+round(redun*omega(eta,b));          
            if k>=1 && k<=neta               
             SST1(k,b)  = SST1(k,b) + STFT(eta,b);
            end
            
%%%%%%%%%%%%%%%%%%%%SST2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
            k = 1+round(redun*omega2(eta,b));
             if k>=1 && k<=neta
               SST2(k,b)  = SST2(k,b) + STFT(eta,b);
             end
             
%%%%%%%%%%%%%%%%%%%%%SST3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
             k = redun*(1+floor(omega3(eta,b)));
             if k>=1 && k<=neta
               SST3(k,b)  = SST3(k,b) + STFT(eta,b);
             end
             
%%%%%%%%%%%%%%%%%%%%%SST4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             k = redun*(1+floor(omega4(eta,b)));
             if k>=1 && k<=neta
               SST4(k,b)  = SST4(k,b) + STFT(eta,b);
             end

        end
    end
end

end

function hp = hermite_poly(t,order)

h{1} = @(x)ones(1,length(x));
h{2} = @(x)2*x;
h{3} = @(x)4*x.^2   - 2;
h{4} = @(x)8*x.^3   - 12*x;
h{5} = @(x)16*x.^4  - 48*x.^2   + 12;
h{6} = @(x)32*x.^5  - 160*x.^3  + 120*x;
h{7} = @(x)64*x.^6  - 480*x.^4  + 720*x.^2      - 120;
h{8} = @(x)128*x.^7 - 1344*x.^5 + 3360*x.^3     - 1680*x;
h{9} = @(x)256*x.^8 - 3584*x.^6 + 13440*x.^4    - 13440*x.^2 + 1680;
h{10}= @(x)512*x.^9 - 9216*x.^7 + 48384*x.^5    - 80640*x.^3 + 30240*x;

hp = h{order+1}(t);

end