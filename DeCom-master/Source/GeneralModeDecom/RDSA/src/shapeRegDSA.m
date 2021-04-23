function [shape,component,SL2,ier] = shapeRegDSA(sig,insAmp,insFreq,insPhase,opt)
% This code implements the algorithm in the paper "A Fast Algorithm for
% Multiresolution Mode Decomposition" by Gao Tang and Haizhao Yang.
%
% Input:
% sig         - given signal
% insAmp    - the instantaneous amplitudes of modes in sig
% insFreq   - the instantaneous frequencies of modes in sig
% insPhase  - the instantaneous phases of modes in sig
% opt       - options for the algorithm
% opt.Ls is the length of the vector representing the shape function
% opt.bandWidth is the band width parameter for the shape function
% opt.diffeoMethod:
%             'intpl', use the spline interpolation to compute
%             diffeomorphisms
%             'nufft', use the NUFFT to compute diffeomorphisms
%
% Output:
% shape     - estimated shape functions
% component - identified modes in sig
% SL2       - the L2-norm of target shape functions
% ier       - if ier<0, the DSA method is not reliable
%           - if ier=0, the DSA method might not be reliable
%
% By Haizhao Yang, 2017

% set up
if ~any(strcmp('diffeoMethod',fieldnames(opt))), opt.diffeoMethod = 'nufft'; end;
if ~any(strcmp('Ls',fieldnames(opt))), opt.Ls = 2000; end;
if ~any(strcmp('bandWidth',fieldnames(opt))), opt.bandWidth = '200'; end;
x_reg = linspace(0,1,opt.Ls);
N = length(sig);
shape = zeros(1,opt.Ls);
component = zeros(1,N);
amp = zeros(1,N);
ier = 1;
uniform_sample = 0:1/N:(1-1/N);
L = round(mean(insFreq));
nonuniform_sample = insPhase/L;

% warpping
switch opt.diffeoMethod
    case 'intpl'% use interpolation
    % 1) there might be an issue at the end points of the interval
    % 2) if the sampling rate is hot high enough, interpolation may not be
    % accurate
    pos1 = find(uniform_sample<nonuniform_sample(1));
    len1 = numel(pos1);
    pos2 = find(uniform_sample>nonuniform_sample(end));
    len2 = numel(pos2);
    if 0
        temp = spline(nonuniform_sample, sig./insAmp, uniform_sample(len1+1:end-len2));
        fpha = sig./insAmp; fpha(len1+1:end-len2) = temp;
    else % periodic extension 
        fpha = spline([uniform_sample(pos1),nonuniform_sample,uniform_sample(pos2)], [sig(end-len1+1:end)./insAmp(end-len1+1:end),sig./insAmp,sig(1:len2)./insAmp(1:len2)], uniform_sample);
    end
    
    % go to the frequency domain to get the Fourier series expansion of the shape
    % function and amplitude function
    fphah = fftshift(fft(fpha))/sqrt(N);
    case 'nufft' % use non-uniform FFT
    % not good if the instantaneous frequencies vary too much, which
    % results in a large condition number for the NUFFT
    fphah=nufft1d1(N,(nonuniform_sample)*2*pi,sig./insAmp,-1,1e-14,N);
    fphah = transpose(fphah)*sqrt(N);
end

% form the matrix fhatMat, which is the matrix H in the RDSA paper
vector1 = (0:L:N/2-L/2)';
vector1 = vector1(1:min(opt.bandWidth+1,numel(vector1)));
vector2 = (-L:-L:-N/2+L/2)';
vector2 = vector2(end:-1:1);
vector2 = vector2(1:min(opt.bandWidth,numel(vector2)));
vector=[vector2;vector1];
K=numel(vector);
u = fphah(vector+floor(N/2)+1);
shapehat =zeros(opt.Ls,1);
if opt.Ls>=K
    st = opt.Ls/2+1-floor(K/2);
    shapehat(st:st+K-1)=u;
else
    st = floor(K/2)+1-floor(opt.Ls/2);
    shapehat=u(st:st+opt.Ls-1);
end
% scaling s.t. shape is independent of opt.Ls
shapeTemp = ifft(fftshift(shapehat))*opt.Ls/sqrt(N);
shape=shape+real(shapeTemp)';
shift = mean(shape);
shape = shape - shift;

if opt.Ls<3*K
    fprintf('Warning: opt.Ls is not large enough,\n');
    fprintf('and hence the RDSA method might not be reliable!\n');
    ier = 0;
end
SL2 = norm(shape)*sqrt(2*pi)/sqrt(length(shape));

% compute the generalized mode
component = insAmp.*spline(x_reg,shape,mod(insPhase,1));
return