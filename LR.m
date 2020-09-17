function [f_est,modos_est,b_e] = LR(f,A,phi,D)
% This function implements the Linear Regression method,
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

b_e = (f*C')/(C*C');

f_est = b_e*C;
modos_est = zeros(I,N);
for i = 1:I
    C = cosenos(A(i,:),phi(i,:),D(i));
    modos_est(i,:) = b_e(2*sum(D(1:i-1))+1:2*sum(D(1:i)))*C;
end;
end