function [shape,component,coef,errorRec,SL2Rec,iter,flag] = shapeDiffusion(sig,numGroup,insAmp,insFreq,insPhase,opt,shapeTrue)
% This code implements the Jacobi style iterative algorithm in "Recursive
% Diffeomorphism-Based Regression for Shape Functions" by Jieren Xu,
% Haizhao Yang, and Ingrid Daubechies, and the Gauss-Seidel style iterative
% scheme proposed in "Multiresolution Mode Decomposition (MMD) for Adaptive Time
% Series Analysis" by Haizhao Yang.
%
% Input:
% sig       - given signal
% numGroup  - the number of generalized mdoes in sig
% insAmp    - the instantaneous amplitudes of modes in sig
% insFreq   - the instantaneous frequencies of modes in sig
% insPhase  - the instantaneous phases of modes in sig
% opt       - options for the algorithm
% opt.maxiter is the maximum iteration number of the iterative scheme
% opt.eps_error is the accuracy parameter for convergence
% opt.eps_diff is the accuracy parameter for convergence
% opt.iterStyle, 'JC' for Jacobi style iteration, 'GS' for Gauss-Seidel
% opt.shapeMethod, method for single shape estimation
%         1: DSA
%         2: BSFK
% opt.show, '1' shows intermediate results, '0' doesn't show
% opt.isCos, '1', compute the hierachical cos expansion coefficient
%            '0', compute the hierachical sin expansion coefficient
% opt.ampFreq, the frequency of the hierachical expansion coefficient
% opt.para contains parameters for estimating a shape function from a
% single mode, say regression, DSA, or other methods.
%
% For the BSFK regression method:
% opt.para.knotremoval_factor is a parameter for the regression method BSFK
% opt.para.nknots is the number of knots for the regression method BSFK
% opt.para.order is the order parameter for the regression method BSFK
% opt.para.Ls is the length of the vector representing the shape function
%
% For the DSA method:
% opt.para.Ls is the length of the vector representing the shape function
% opt.para.bandWidth is the band width parameter for the shape function
% opt.para.diffeoMethod:
%             'intpl', use the spline interpolation to compute
%             diffeomorphisms
%             'nufft', use the NUFFT to compute diffeomorphisms
%
% Optional inputs:
% shapeTrue - the ground true shapes for the purpose of debuging
%
% Output:
% shape     - estimated shape functions
% component - identified modes in sig
% coef      - multiresolution expansion coefficients
% errorRec  - a record of errors in each iteration
% SL2Rec    - a record of the L2-norm of target shape functions in each
%             iteration
% iter      - the number of iterations needed for convergence
% flag      - a variable containing the information about the convergence
%
% By Haizhao Yang, 2017


if ~any(strcmp('maxiter',fieldnames(opt))), opt.maxiter = 200; end;
if ~any(strcmp('eps_error',fieldnames(opt))), opt.eps_error = 1e-6; end;
if ~any(strcmp('show',fieldnames(opt))), opt.show = 0; end;
if ~any(strcmp('iterStyle',fieldnames(opt))), opt.iterStyle = 'GS'; end;
if ~any(strcmp('shapeMethod',fieldnames(opt))), opt.shapeMethod = 1; end;
if ~any(strcmp('eps_diff',fieldnames(opt))), opt.eps_diff = 1e-6; end;
if ~any(strcmp('para',fieldnames(opt))), opt.para = []; end;
if ~any(strcmp('isCos',fieldnames(opt))), opt.isCos = 1; end;
if ~any(strcmp('ampFreq',fieldnames(opt))), opt.ampFreq = 0; end;

%% set up
maxiter = opt.maxiter;
eps_error = opt.eps_error;
shape = cell(1,numGroup);
coef = cell(1,numGroup);
prefac = cell(1,numGroup);
shapeInLoop = cell(1,numGroup);
componentInLoop = cell(1,numGroup);
component = cell(1,numGroup);

%% initialization
iter = 0;
SL2Rec = [];
isBreak = 0;
delta = 1;
error = norm(sig)/sqrt(length(sig));
errorInit = error;
errorRec = error;
errorOld = 2*error;
r_n = sig;
for cntGroup = 1:numGroup
    shape{cntGroup} = zeros(1,opt.para.Ls);
    component{cntGroup} = zeros(1,length(sig));
    L = round(mean(insFreq(cntGroup,:)));
    if opt.isCos == 1
        prefac{cntGroup} = cos(2*pi*opt.ampFreq*mod(insPhase(cntGroup,:)/L,1));
    else
        prefac{cntGroup} = sin(2*pi*opt.ampFreq*mod(insPhase(cntGroup,:)/L,1));
    end
end

if opt.shapeMethod == 1
    if numel(opt.para.bandWidth) > 1
        bandWidthbk = opt.para.bandWidth;
    else
        bandWidthbk = [];
    end
end

%% self-consistent iteration
while iter< maxiter && error> eps_error && delta>eps_error && ...
        abs(errorOld-error)>opt.eps_diff
    % in the case of clean signal and bad parameters of regression, this can stop redundant iterations
    iter = iter + 1;
    f = r_n;
    errorOld = error;
    
    SL2 = 0;
    delta = 0;
    for cntGroup = 1:numGroup
        switch opt.shapeMethod
            case 1
                if numel(bandWidthbk) > 0
                    opt.para.bandWidth = bandWidthbk(cntGroup);
                end
                [shapeInLoop{cntGroup},componentInLoop{cntGroup},SL2temp,ier] = shapeRegDSA(prefac{cntGroup}.*f,insAmp(cntGroup,:),insFreq(cntGroup,:),insPhase(cntGroup,:),opt.para);
            case 2
                [shapeInLoop{cntGroup},componentInLoop{cntGroup},SL2temp,ier] = shapeRegBSFK(prefac{cntGroup}.*f,insAmp(cntGroup,:),insPhase(cntGroup,:),opt.para);
        end
        if opt.ampFreq ~=0
            componentInLoop{cntGroup} = prefac{cntGroup}.*componentInLoop{cntGroup}*2;
            shapeInLoop{cntGroup} = shapeInLoop{cntGroup}*2;
        end
        SL2 = max(SL2,SL2temp);
        if ier<0
            isBreak = 1;
        end
        shape{cntGroup} = shape{cntGroup} + shapeInLoop{cntGroup};
        component{cntGroup} = component{cntGroup} + componentInLoop{cntGroup};
        delta = max(delta,norm(componentInLoop{cntGroup})/sqrt(length(componentInLoop{cntGroup})));
        r_n  = r_n - componentInLoop{cntGroup};
        switch opt.iterStyle
            case 'JC'
                % do nothing
            case 'GS'
                f = r_n;
        end
    end
    SL2Rec = [SL2Rec SL2];
    %% compute error
    error = norm(r_n)/sqrt(length(r_n));
    if ier <= 0 && errorOld<error  % safegard
        isBreak = 1;
    end
    errorRec = [errorRec error];
    if error>errorInit % safegard
        isBreak = 1;
    end
    %iter
    %[SL2 delta error]
    if isBreak
        break;
    end
end

for cntGroup = 1:numGroup
    coef{cntGroup} = sqrt(sum(abs(shape{cntGroup}).^2)*sqrt(2*pi)/length(shape{cntGroup}));
    shape{cntGroup} =  shape{cntGroup}/coef{cntGroup};
end

if opt.show  && nargin > 5
    x_reg = linspace(0,1,opt.para.Ls);
    scrsz = get(groot,'ScreenSize');
    for i = 1:numGroup
        fig_opts{i}=struct('Position',[scrsz(4)/2   scrsz(4)-(scrsz(4)/2)*i    scrsz(4)/2   scrsz(4)/3] );
        fig_opts_Initial{i} = struct('Position',[scrsz(4)   scrsz(4)-(scrsz(4)/2)*i    scrsz(4)/2   scrsz(4)/3] );
        Fig_opts{i}=struct('Position',[20   scrsz(4)-(scrsz(4)/2)*i    scrsz(4)/2   scrsz(4)/3] );
    end
    for i = 1: numGroup
        fig = figure(i+10);set(fig,fig_opts{i});
        subplot(211);
        plot(x_reg,shape{i},'linewidth',2); hold on;
        plot(x_reg,shapeTrue{i}(x_reg),'linewidth',2);hold off;
        legend('reg','true')
        title(['iter',num2str(iter),' for s',num2str(i)])
        subplot(212);
    end
end
flag = zeros(1,5);
if abs(errorOld-error)<=opt.eps_diff
    if opt.show
        fprintf('abs(diff)<tol\n');
    end
    flag(1) = 1;
end
if error>errorOld
    if opt.show
        fprintf('residual norm grows up\n');
    end
    flag(2) = 1;
end
if error<=eps_error | delta<=eps_error
    if opt.show
        fprintf('residual norm is smaller than tolerance\n');
    end
    flag(3) = 1;
end
if iter>=maxiter
    if opt.show
        fprintf('number of iterations > maxiter\n');
    end
    flag(4) = 1;
end
if isBreak;
    if opt.show
        fprintf('shape function extraction is not reliable\n');
    end
    flag(5) = 1;
end




