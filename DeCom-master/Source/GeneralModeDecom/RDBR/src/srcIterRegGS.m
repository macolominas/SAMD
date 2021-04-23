function [shape,component,errorRec,SL2Rec,iter,flag] = srcIterRegGS(fff,N,numGroup,insAmp,insPhase,opt,fTrue,shapeTrue)
% This code implements the iterative regression algorithm in "Recursive 
% Diffeomorphism-Based Regression for Shape Functions" by Jieren Xu,
% Haizhao Yang, and Ingrid Daubechies. We adopt the Gauss-Seidel iterative 
% scheme proposed in "Multiresolution Mode Decomposition for Adaptive Time 
% Series Analysis" by Haizhao Yang.
%
% Input:
% fff       - given signal
% N         - length of fff
% numGroup  - the number of generalized mdoes in fff
% insAmp    - the instantaneous amplitudes of modes in fff
% insPhase  - the instantaneous phases of modes in fff
% opt       - options for the algorithm
% opt.maxiter is the maximum iteration number of the iterative scheme
% opt.eps_error is the accuracy parameter for convergence
% opt.eps_diff is the accuracy parameter for convergence
% opt.knotremoval_factor is a parameter for the regression method BSFK
% opt.nknots is the number of knots for the regression method BSFK
% opt.order is the order parameter for the regression method BSFK
% opt.show, '1' shows intermediate results, '0' doesn't show
%
% Optional inputs:
% fTrue     - the ground true modes in fff for the purpose of debuging
% shapeTrue - the ground true shapes for the purpose of debuging
%
% Output:
% shape     - estimated shape functions
% component - identified modes in fff
% errorRec  - a record of errors in each iteration
% SL2Rec    - a record of the L2-norm of target shape functions in each
%             iteration
% iter      - the number of iterations needed for convergence
% flag      - a variable containing the information about the convergence
%
% By Haizhao Yang

%% set up
maxiter = opt.maxiter;
eps_error = opt.eps_error;
options = struct('animation', 0, 'knotremoval_factor',opt.knotremoval_factor);
if opt.show
    scrsz = get(groot,'ScreenSize');
    for i = 1:numGroup
        fig_opts{i}=struct('Position',[scrsz(4)/2   scrsz(4)-(scrsz(4)/2)*i    scrsz(4)/2   scrsz(4)/3] );
        fig_opts_Initial{i} = struct('Position',[scrsz(4)   scrsz(4)-(scrsz(4)/2)*i    scrsz(4)/2   scrsz(4)/3] );
        Fig_opts{i}=struct('Position',[20   scrsz(4)-(scrsz(4)/2)*i    scrsz(4)/2   scrsz(4)/3] );
    end
end
pp = cell(1,numGroup);
shape = cell(1,numGroup);
shapeInLoop = cell(1,numGroup);
componentInLoop = cell(1,numGroup);
component = cell(1,numGroup);

%% initialization
nknots = opt.nknots;
order = opt.order;
x_reg = linspace(0,1,1000)';
std = 0;
iter = 0;
r_n = fff;
SL2 = 0;
for cntGroup = 1:numGroup
    phi2 = mod(insPhase(cntGroup,:), 1);
    X = phi2(:);   Y = r_n./insAmp(cntGroup,:); Y = Y(:);
    [pp{cntGroup}] = BSFK(X',Y,order,nknots,[],options);
    Yest = ppval(pp{cntGroup},X);
    std = max(std,norm(Y-Yest)/sqrt(length(Y)));

    shape{cntGroup} = ppval(pp{cntGroup},x_reg);
    SL2 = max(SL2,norm(shape{cntGroup})*sqrt(2*pi)/sqrt(length(shape{cntGroup})));
    shift = mean(shape{cntGroup});
    shape{cntGroup} = shape{cntGroup} - shift;
    componentInLoop{cntGroup} = (ppval(pp{cntGroup},mod(insPhase(cntGroup,:),1))-shift).*insAmp(cntGroup,:);
    r_n  = r_n - componentInLoop{cntGroup};
    component{cntGroup} = componentInLoop{cntGroup};
end
SL2Rec = SL2;


%% self-consistent iteration
error = norm(r_n)/sqrt(length(r_n));
errorRec = error;
old_error = 2*error;
% display iteration
disp(['initialization error is  ',num2str(error)])


if opt.show && nargin > 6
    show_ind = 1:min(N,200);
    for i = 1:numGroup
        fig = figure(i+1000);set(fig,fig_opts_Initial{i});   % initial estimation
        % show shapeInLoop
        subplot(211);
        plot(x_reg,ppval(pp{i},x_reg),'linewidth',2); hold on;
        plot(x_reg,shapeTrue{i}(x_reg),'linewidth',2);hold off;
        legend('reg','true')
        title(['initial for s',num2str(i)])
        subplot(212);
        % plot(componentInLoop{i}(show_ind));hold on;plot(fTrue{i}(show_ind));hold off;legend('initial','true')
    end
end

isBreak = 0;
delta = 1;
while iter< maxiter && error> eps_error && delta>eps_error && ...
        abs(old_error-error)>opt.eps_diff
    % in the case of clean signal and bad parameters of regression, this can stop redundant iterations
    iter = iter + 1;
   % f = r_n;
    old_error = error;
    
    % extract s1 and s2 from f
    SL2 = 0;
    std = 0;
    delta = 0;
    for cntGroup = 1:numGroup
        % warping and folding
        phi2 = mod(insPhase(cntGroup,:), 1) ;
        X = phi2(:);    Y = r_n./insAmp(cntGroup,:); Y = Y(:);
        
        % BSFK
        [pp{cntGroup},ier] = BSFK(X',Y,order,nknots,[],options);
        if ier<0
            isBreak = 1;
        end
        
        % get shape
        shapeInLoop{cntGroup} = ppval(pp{cntGroup},x_reg);
        Yest = ppval(pp{cntGroup},X);
        std = max(std,norm(Y-Yest)/sqrt(length(Y)));
        SL2 = max(SL2,norm(shapeInLoop{cntGroup})*sqrt(2*pi)/sqrt(length(shapeInLoop{cntGroup})));
        shift = mean(shapeInLoop{cntGroup}); % this can avoid breaking down
        shapeInLoop{cntGroup} = shapeInLoop{cntGroup} - shift;
        
        % update iteration
        shape{cntGroup} = shape{cntGroup} + shapeInLoop{cntGroup};
        componentInLoop{cntGroup} = (ppval(pp{cntGroup},mod(insPhase(cntGroup,:),1))-shift).*insAmp(cntGroup,:);
        delta = max(delta,norm(componentInLoop{cntGroup})/sqrt(length(componentInLoop{cntGroup})));
        component{cntGroup} = component{cntGroup} + componentInLoop{cntGroup};
        r_n  = r_n - componentInLoop{cntGroup};
    end
    SL2Rec = [SL2Rec SL2];
    %% compute error
    error = norm(r_n)/sqrt(length(r_n));
    errorRec = [errorRec error];
    iter
    [std^2 SL2 delta error]
    if isBreak %| SL2>SL2old
        break;
    end
end

if opt.show  && nargin > 6
    for i = 1: numGroup
        fig = figure(i+10);set(fig,fig_opts{i});
        % show shapeInLoop
        subplot(211);
        plot(x_reg,shape{i},'linewidth',2); hold on;
        plot(x_reg,shapeTrue{i}(x_reg),'linewidth',2);hold off;
        legend('reg','true')
        title(['iter',num2str(iter),' for s',num2str(i)])
        subplot(212);
        %   plot(component{i}( show_ind ));hold on;plot(fTrue{i}( show_ind ));hold off; legend('iter','true')
    end
end
flag = zeros(1,5);
if abs(old_error-error)<opt.eps_diff; fprintf('abs(diff)<tol\n'); flag(1) = 1; end;
if error>old_error; fprintf('residual norm grows up\n'); flag(2) = 1; end;
if error<eps_error; fprintf('residual norm is smaller than tolerance\n'); flag(3) = 1; end;
if iter>=maxiter;  fprintf('number of iterations > maxiter\n'); flag(4) = 1; end;
if isBreak; fprintf('regression is not reliable\n'); flag(5) = 1; end;




