function [f_est,b_e] = one_wsf_nonlinear(f,A1,phi1,D1,order_Phi)

C = [cosenos(A1,phi1,D1)];
coef_est = (f*C')/(C*C');

% A1 = X(1,:);
% phi1 = X(2,:);
% A2 = X(3,:);
% phi2 = X(4,:);

% c1 = b(1:D1);
% d1 = b(D1+1:2*D1);
% alpha1 = b(2*D1+1:end);
% 
% alpha1_inic = zeros(order_Phi,D1);
% alpha1_inic(1,:) = 1:D1;

modelfun = @(b,X)(...
    X(1,:).*(b(1:D1)*cos(reshape([1 zeros(1,order_Phi-1) b(2*D1+1:end)],[order_Phi,D1])'*phi_aux(X(2,:),order_Phi)) ...
    +   b(D1+1:2*D1)*sin(reshape([1 zeros(1,order_Phi-1) b(2*D1+1:end)],[order_Phi,D1])'*phi_aux(X(2,:),order_Phi))));

alpha1_inic = zeros(order_Phi,D1-1);
alpha1_inic(1,:) = 2:D1;

valor_inicial = [coef_est alpha1_inic(:)'];

% [b_e,R,J,CovB,MSE,ErrorModelInfo] = nlinfit([A1;phi1;A2;phi2],f,modelfun,valor_inicial)
options = statset;
%options.RobustWgtFun = 'cauchy'; %sin esto no anda en muchos ejemplos. En
% el ejemplo de tren de gaussianas delgadas no anda. El algoritmo no itera:
% toma como solución el valor inicial.
options.RobustWgtFun = 'cauchy';
options.MaxIter = 3000;
% options.TolX = 1e-10;
b_e = nlinfit([A1;phi1],f,modelfun,valor_inicial,options);

f_est = modelfun(b_e,[A1;phi1]);
end

function M = phi_aux(phi,order)
  largo = length(phi);
  M = zeros(order,largo);
  for ii = 1:order
      M(ii,:) = phi.^(ii)/(factorial(ii));
  end;
end