function [f_est,modos_est,b_e] = SAMD(f,A,phi,D,order_Phi)
% This function implements the Shape-adaptive mode decomposition method,
% from the paper "Signal decomposition for time-varying wave-shape
% functions and its biomedical application" by Marcelo A. Colominas and
% Hau-Tieng Wu.
%
% SYNTAX
%
% INPUT:
% f: row vector containing the signal to decompose
% A: matrix containing the instantanaeous amplitudes (one per row)
% phi: matrix containing the instantanaeous phases (one per row)
% D: vector containing the parameters D_i.
% order_Phi: scalar defining the order of estimation for the phases.
%
% OUTPUT:
% f_est: estimated reconstructed signal
% modos_est: estimated components (one per row) (same size as A and phi)
% b_e: regression coefficients vector
%
% Marcelo A. Colominas
% email: macolominas@conicet.gov.ar
% 17-SEP-2020

[I, N] = size(A);
C = zeros(2*sum(D),N);
for i = 1:I
    C(2*sum(D(1:i-1))+1:2*sum(D(1:i)),:) = cosenos(A(i,:),phi(i,:),D(i));
end

coef_est = (f*C')/(C*C');

modelfun = @(b,X)(regresion(b,X,D,order_Phi));

valor_inicial = coef_est;

for i = 1:I
    alpha_inic = zeros(order_Phi,D(i)-1);
    alpha_inic(1,:) = 2:D(i);
    valor_inicial = [valor_inicial alpha_inic(:)'];
end

options = statset;
options.RobustWgtFun = 'cauchy';
options.MaxIter = 3000;

warning('off')
b_e = nlinfit([A;phi],f,modelfun,valor_inicial,options);
warning('on')

f_est = modelfun(b_e,[A;phi]);
modos_est = zeros(I,N);
for i = 1:I
    A_aux = zeros(I,N);
    A_aux(i,:) = A(i,:);
    phi_ = zeros(I,N);
    phi_(i,:) = phi(i,:);
    modos_est(i,:) = modelfun(b_e,[A_aux;phi_]);
end
end

function v = regresion(b,X,D,order_Phi)
    I = length(D);
    A = X(1:I,:);
    phi = X(I+1:end,:);
    
    v = A(1,:).*(b(1:D(1))*cos(reshape([1 zeros(1,order_Phi-1) b(2*(sum(D(1:I)))+1:2*(sum(D(1:I)))+(D(1)-1)*order_Phi)],[order_Phi,D(1)])'*phi_aux(phi(1,:),order_Phi)) ...
    +   b(D(1)+1:2*D(1))*sin(reshape([1 zeros(1,order_Phi-1) b(2*(sum(D(1:I)))+1:2*(sum(D(1:I)))+(D(1)-1)*order_Phi)],[order_Phi,D(1)])'*phi_aux(phi(1,:),order_Phi)));
    
    for i = 2:I
        v = v + A(i,:).*(b(2*sum(D(1:i-1))+1:2*sum(D(1:i-1))+D(i))*cos(reshape([1 zeros(1,order_Phi-1) b(2*sum(D(1:I))+(sum(D(1:i-1))-i+1)*order_Phi+1:2*sum(D(1:I))+(sum(D(1:i))-i)*order_Phi)],[order_Phi,D(i)])'*phi_aux(phi(i,:),order_Phi)) ...
                     +     b(2*sum(D(1:i-1))+D(i)+1:2*sum(D(1:i)))*sin(reshape([1 zeros(1,order_Phi-1) b(2*sum(D(1:I))+(sum(D(1:i-1))-i+1)*order_Phi+1:2*sum(D(1:I))+(sum(D(1:i))-i)*order_Phi)],[order_Phi,D(i)])'*phi_aux(phi(i,:),order_Phi)));
    end  
end

function M = phi_aux(phi,order)
  largo = length(phi);
  M = zeros(order,largo);
  for ii = 1:order
      M(ii,:) = phi.^(ii)/(factorial(ii));
  end;
end