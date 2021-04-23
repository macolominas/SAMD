function [f_est,WSF1,coef_est] = one_wsf_linear(f,A1,phi1,D1)
C = cosenos(A1,phi1,D1);
C1_aux = cosenos(ones(1,1000),-pi:2*pi/1000:pi-2*pi/1000,D1);

coef_est = (f*C')/(C*C');
f_est = coef_est*C;

WSF1 = coef_est*C1_aux;