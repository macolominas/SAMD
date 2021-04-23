function [f_est,f1_est,f2_est,WSF1,WSF2,coef_est] = two_wsf_linear(f,A1,phi1,D1,A2,phi2,D2)
C = [cosenos(A1,phi1,D1);cosenos(A2,phi2,D2)];
C1_aux = cosenos(ones(1,1000),-pi:2*pi/1000:pi-2*pi/1000,D1);
C2_aux = cosenos(ones(1,1000),-pi:2*pi/1000:pi-2*pi/1000,D2);

coef_est = (f*C')/(C*C');
f_est = coef_est*C;
f1_est = coef_est(1:2*D1)*C(1:2*D1,:);
f2_est = coef_est(2*D1+1:end)*C(2*D1+1:end,:);

WSF1 = coef_est(1:2*D1)*C1_aux;
WSF2 = coef_est(2*D1+1:end)*C2_aux;  