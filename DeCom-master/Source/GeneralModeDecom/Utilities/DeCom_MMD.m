function [shape,component,Hcoef,flag,idx,iter,iterDR] = DeCom_MMD(sig,x,numGroup,insAmp,insFreq,insPhase,opt)
% This code implements the robust iterative regression algorithm in
% "Multiresolution Mode Decomposition (MMD) for Adaptive Time Series Analysis"
% by Haizhao Yang.
%
% Input:
% sig       - a given real signal
% x         - the time position of sig, i.e. sig = f(x) for some function f
% numGroup  - the number of generalized mdoes in sig
% insAmp    - the instantaneous amplitudes of modes in sig, a matrix with
%             entires all 1
% insFreq   - the instantaneous frequencies of modes in sig
% insPhase  - the instantaneous phases of modes in sig
% opt       - options for the algorithm
% opt.numSweep is the maximum number of sweeping in the multiresolution
% frequency domain
% opt.accuracy is the accuracy parameter for the convergence of MMD
% opt.eps_diff is the accuracy parameter for convergence
% opt.ampErrBandWidth is the band width of the amplitude estimation error in
% the frequency domain.
% opt.shapeMethod, method for single shape estimation
%         1: DSA
%         2: BSFK
% opt.lowestCompEng is the lowest energy for a component, i.e. the L2-norm
% of a component is larger than or equal to opt.lowestCompEng.
% opt.decayBndExpCoef is the bound of the decay rate of the multiresolution
% expansion coefficients, i.e., |a_{n+1}|/|a_{n}| <=  opt.decayBndExpCoef
% and |b_{n+1}|/|b_{n}| <=  opt.decayBndExpCoef.
% opt.para contains parameters for estimating a shape function from a
% single mode, say regression, DSA, or other methods.
%
%
% opt.para.show, '1' shows intermediate results, '0' doesn't show
% opt.para.maxiter is the maximum iteration number of the iterative scheme
% opt.para.eps_error is the accuracy parameter for convergence
% opt.para.eps_diff is the accuracy parameter for convergence
% opt.para.iterStyle, 'JC' for Jacobi style iteration, 'GS' for Gauss-Seidel
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
% Output:
% shape     - estimated shape functions
% component - identified modes in sig
% Hcoef     - a cell variable with Hcoef{i} as a structure storing the
%             multiresolution expansion coefficients for the ith component;
%             Hcoef{i}.a0 for a_{0}, Hcoef{i}.a(n) for a_{n} and
%             Hcoef{i}.b(n) for b_{n}.
% flag(1)   - 0, MMD achieves target depth; 1, MMD stops with a smaller
%             depth
% flag(2)   - 0, numGroup is equal to the exact number of components
%             -1, numGroup is smaller than the exact number of components
%             n, n>0, numGroup is larger than the exact number of
%             components; the new estimation of components is n.
% idx       - the indices of significant components above
%             opt.lowestCompEng.
% iter      - the number of sweeping in practice
% iterDR - the total number of diffeomorphism regression
%
% By Haizhao Yang, 2017

if ~any(strcmp('ampErrBandWidth',fieldnames(opt))), opt.ampErrBandWidth = '20'; end;
if ~any(strcmp('lowestCompEng',fieldnames(opt))), opt.lowestCompEng = 0; end;
if ~any(strcmp('decayBndExpCoef',fieldnames(opt))), opt.decayBndExpCoef = 1e10; end;
if ~any(strcmp('numSweep',fieldnames(opt))), opt.numSweep = 10; end;
if ~any(strcmp('accuracy',fieldnames(opt))), opt.accuracy = 1e-6; end;
if ~any(strcmp('eps_diff',fieldnames(opt))), opt.eps_diff = 1e-6; end;

% initialization
componentIn = cell(1,numGroup);
errorOld = 1e10;
errSweepOld = 1e10;
flag = [0,0];
isCoefBreak = 0;
iterDR = 0;
N = numel(sig);
Hcoef = cell(1,numGroup);
for cntc = 1:numGroup
    Hcoef{cntc}.a = zeros(1,2*opt.ampErrBandWidth);
    Hcoef{cntc}.b = zeros(1,2*opt.ampErrBandWidth);
end
shape = cell(1,numGroup);
for cntc = 1:numGroup
    shape{cntc}.scn = cell(1,2*opt.ampErrBandWidth);
    shape{cntc}.ssn = cell(1,2*opt.ampErrBandWidth);
    shape{cntc}.s0 = zeros(1,opt.para.Ls);
    for cnt = 1:2*opt.ampErrBandWidth
        shape{cntc}.scn{cnt} = zeros(1,opt.para.Ls);
        shape{cntc}.ssn{cnt} = zeros(1,opt.para.Ls);
    end
end
idx = 1:numGroup;
vec = 1:opt.ampErrBandWidth;
if opt.ampErrBandWidth > 0
    vec = [vec;-vec];
    vec = reshape(vec,[1,2*size(vec,2)]);
end
vec = [0,vec];
numSweep = opt.numSweep;
normSig = norm(sig);
for cnts = 1:numSweep
    for cnt = 1:numel(vec)
        opt.isCos = 1; opt.ampFreq = vec(cnt);
        [shapeIn,compTemp,coef,errorRec,~,iter] = shapeDiffusion(sig,numGroup,insAmp,insFreq,insPhase,opt);
        iterDR = iterDR + iter*numGroup;
        if vec(cnt) == 0
            error = errorRec(end);%todo
        else
            errorCos = errorRec(end);
        end
        for cntc = 1:numGroup
            if vec(cnt) == 0
                componentIn{cntc} = compTemp{cntc};
                shape{cntc}.s0 = shape{cntc}.s0 + shapeIn{cntc}*coef{cntc};
            else
                componentIn{cntc} = componentIn{cntc} + compTemp{cntc};
                shape{cntc}.scn{cnt-1} = shape{cntc}.scn{cnt-1} + shapeIn{cntc}*coef{cntc};
                
                if vec(cnt) == 1
                    if norm(shape{cntc}.scn{cnt-1})>norm(shape{cntc}.s0)*opt.decayBndExpCoef
                        isCoefBreak = 1;
                        break;
                    end
                end
            end
            sig = sig - compTemp{cntc};
        end
        if isCoefBreak
            break;
        end
        if vec(cnt)~=0
            opt.isCos = 0; opt.ampFreq = vec(cnt);
            [shapeIn,compTemp,coef,errorRec,~,iter] = shapeDiffusion(sig,numGroup,insAmp,insFreq,insPhase,opt);
            iterDR = iterDR + iter*numGroup;
            for cntc = 1:numGroup;
                componentIn{cntc} = componentIn{cntc} + compTemp{cntc};
                shape{cntc}.ssn{cnt-1} = shape{cntc}.ssn{cnt-1} + shapeIn{cntc}*coef{cntc};
                sig = sig - compTemp{cntc};
                
                if vec(cnt) == 1
                    if norm(shape{cntc}.ssn{cnt-1})>norm(shape{cntc}.s0)*opt.decayBndExpCoef
                        isCoefBreak = 1;
                        break;
                    end
                end
            end
            error = (errorCos+errorRec(end))/2;
        end
        if isCoefBreak
            break;
        end
        if errorOld < error
            flag(1) = 1;
            break;
        else
            errorOld = error;
        end
    end
    if cnts == 1
        component = componentIn;
    else
        for cntg = 1:numGroup
            component{cntg} = component{cntg} + componentIn{cntg};
        end
    end
    cnts
    errSweep = norm(sig)/normSig
    if errSweep > errSweepOld - opt.eps_diff, break; end; % cannot improve any more
    errSweepOld = errSweep;
    if errSweep < opt.accuracy, break; end;
end
iter = cnts;

% compute the multiresolution coefficients
for cnt = 1:numel(vec)
    for cntc = 1:numGroup
        if vec(cnt) == 0
            Hcoef{cntc}.a0 = sqrt(sum(abs(shape{cntc}.s0).^2)/length(shape{cntc}.s0));
            shape{cntc}.s0 = shape{cntc}.s0/Hcoef{cntc}.a0;
        else
            Hcoef{cntc}.a(cnt-1) = sqrt(sum(abs(shape{cntc}.scn{cnt-1}).^2)/length(shape{cntc}.scn{cnt-1}));
            shape{cntc}.scn{cnt-1} = shape{cntc}.scn{cnt-1}/Hcoef{cntc}.a(cnt-1);
        end
    end
    if vec(cnt)~=0
        for cntc = 1:numGroup;
            Hcoef{cntc}.b(cnt-1) = sqrt(sum(abs(shape{cntc}.ssn{cnt-1}).^2)/length(shape{cntc}.ssn{cnt-1}));
            shape{cntc}.ssn{cnt-1} = shape{cntc}.ssn{cnt-1}/Hcoef{cntc}.b(cnt-1);
        end
    end
end

if isCoefBreak > 0
    flag(2) = -1;
    idx = 1:numGroup;
end
