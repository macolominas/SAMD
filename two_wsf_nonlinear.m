function [f_est,f1_est,f2_est,b_e] = two_wsf_nonlinear(f,A1,phi1,D1,A2,phi2,D2,order_Phi)

C = [cosenos(A1,phi1,D1);cosenos(A2,phi2,D2)];
coef_est = (f*C')/(C*C');

% A1 = X(1,:);
% phi1 = X(2,:);
% A2 = X(3,:);
% phi2 = X(4,:);

% c1 = b(1:D1);
% d1 = b(D1+1:2*D1);
% c2 = b(2*D1+1:2*D1+D2);
% d2 = b(2*D1+D2+1:2*(D1+D2));
% alpha1 = b(2*(D1+D2)+1:2*(D1+D2)+D1*order_Phi);
% alpha2 = b(2*(D1+D2)+D1*order_Phi+1:end);

% modelfun = @(b,X)(...
%     X(1,:).*(b(1:D1)*cos(reshape(b(2*(D1+D2)+1:2*(D1+D2)+D1*order_Phi),[order_Phi,D1])'*phi_aux(X(2,:),order_Phi)) ...
%     +   b(D1+1:2*D1)*sin(reshape(b(2*(D1+D2)+1:2*(D1+D2)+D1*order_Phi),[order_Phi,D1])'*phi_aux(X(2,:),order_Phi)))...
%     + X(3,:).*(b(2*D1+1:2*D1+D2)*cos(reshape(b(2*(D1+D2)+D1*order_Phi+1:end),[order_Phi,D2])'*phi_aux(X(4,:),order_Phi)) ...
%     +     b(2*D1+D2+1:2*(D1+D2))*sin(reshape(b(2*(D1+D2)+D1*order_Phi+1:end),[order_Phi,D2])'*phi_aux(X(4,:),order_Phi))));
% 
% alpha1_inic = zeros(order_Phi,D1);
% alpha1_inic(1,:) = 1:D1;
% 
% alpha2_inic = zeros(order_Phi,D2);
% alpha2_inic(1,:) = 1:D2;

modelfun = @(b,X)(...
    X(1,:).*(b(1:D1)*cos(reshape([1 zeros(1,order_Phi-1) b(2*(D1+D2)+1:2*(D1+D2)+(D1-1)*order_Phi)],[order_Phi,D1])'*phi_aux(X(2,:),order_Phi)) ...
    +   b(D1+1:2*D1)*sin(reshape([1 zeros(1,order_Phi-1) b(2*(D1+D2)+1:2*(D1+D2)+(D1-1)*order_Phi)],[order_Phi,D1])'*phi_aux(X(2,:),order_Phi)))...
    + X(3,:).*(b(2*D1+1:2*D1+D2)*cos(reshape([1 zeros(1,order_Phi-1) b(2*(D1+D2)+(D1-1)*order_Phi+1:end)],[order_Phi,D2])'*phi_aux(X(4,:),order_Phi)) ...
    +   b(2*D1+D2+1:2*(D1+D2))*sin(reshape([1 zeros(1,order_Phi-1) b(2*(D1+D2)+(D1-1)*order_Phi+1:end)],[order_Phi,D2])'*phi_aux(X(4,:),order_Phi))));

alpha1_inic = zeros(order_Phi,D1-1);
alpha1_inic(1,:) = 2:D1;

alpha2_inic = zeros(order_Phi,D2-1);
alpha2_inic(1,:) = 2:D2;

valor_inicial = [coef_est alpha1_inic(:)' alpha2_inic(:)'];
% valor_inicial = [1 .5 0 0 .75 1 0 0 1.98 0 0 2.05 0.0002 0]

% [b_e,R,J,CovB,MSE,ErrorModelInfo] = nlinfit([A1;phi1;A2;phi2],f,modelfun,valor_inicial)
options = statset;
%options.RobustWgtFun = 'cauchy'; %sin esto no anda en muchos ejemplos. En
% el ejemplo de tren de gaussianas delgadas no anda. El algoritmo no itera:
% toma como solución el valor inicial.
options.RobustWgtFun = 'cauchy';
options.MaxIter = 3000;
% options.TolX = 1e-10;
b_e = nlinfit([A1;phi1;A2;phi2],f,modelfun,valor_inicial,options);

f_est = modelfun(b_e,[A1;phi1;A2;phi2]);
f1_est = modelfun(b_e,[A1;phi1;zeros(size(A1));zeros(size(A1))]);
f2_est = modelfun(b_e,[zeros(size(A2));zeros(size(A2));A2;phi2]);
end

function M = phi_aux(phi,order)
  largo = length(phi);
  M = zeros(order,largo);
  for ii = 1:order
      M(ii,:) = phi.^(ii)/(factorial(ii));
  end;
end