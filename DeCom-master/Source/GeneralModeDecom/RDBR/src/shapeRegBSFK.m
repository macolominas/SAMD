function [shape,component,SL2,ier] = shapeRegBSFK(sig,insAmplitude,insPhase,opt)
% This code uses the BSFK method to estimate the shape function of a single
% generalized mode.
%
% Input:
% sig         - given signal
% insAmp    - the instantaneous amplitudes of modes in sig
% insPhase  - the instantaneous phases of modes in sig
% opt       - options for the algorithm
% opt.knotremoval_factor is a parameter for the regression method BSFK
% opt.nknots is the number of knots for the regression method BSFK
% opt.order is the order parameter for the regression method BSFK
% opt.Ls is the length of the vector representing the shape function
%
% Output:
% shape     - estimated shape functions
% component - identified modes in sig
% SL2       - the L2-norm of target shape functions
% ier       - if ier<0, regression is not reliable
%
% By Haizhao Yang, 2017


if ~any(strcmp('knotremoval_factor',fieldnames(opt))), opt.knotremoval_factor = 1.01; end;
if ~any(strcmp('nknots',fieldnames(opt))), opt.nknots = 20; end;
if ~any(strcmp('order',fieldnames(opt))), opt.order = 3; end;
if ~any(strcmp('Ls',fieldnames(opt))), opt.Ls = 1000; end;

%% set up
options = struct('animation', 0, 'knotremoval_factor',opt.knotremoval_factor);
nknots = opt.nknots;
order = opt.order;

x_reg = linspace(0,1,opt.Ls);
phi2 = mod(insPhase, 1);
X = phi2(:);   Y = sig./insAmplitude; Y = Y(:);
[pp,ier] = BSFK(X',Y,order,nknots,[],options);
shape = ppval(pp,x_reg);
shift = mean(shape);
shape = shape - shift;
SL2 = norm(shape)*sqrt(2*pi)/sqrt(length(shape));
component = (ppval(pp,mod(insPhase,1))-shift).*insAmplitude;