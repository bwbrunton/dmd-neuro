function [Xhat z0] = DMD_recon(Phi, lambda, x0, M, varargin)
% function [Xhat z0] = DMD_recon(Phi, lambda, x0, M, varargin)
% reconstructs the data matrix X from its Dynamic Mode Decomposition,
% where [Phi mu lambda diagS x0] = DMD(X)
%
% INPUTS: 
%
% Phi, lambda, x0 are outputs from DMD.m
% M is the number of time points desired in the reconstruction
%
% optional parameters: {'parameter_name', [default_value]}
%   {'keep_modes', []}, if empty, use all modes in reconstruction of Xhat
%                       if not empty, use only modes indexed here in the
%                       reconstruction of Xhat
%
%
% OUTPUTS:
%
% Xhat, the reconstructed data matrix, will be n x M where
%       n is the length of x0
%
%       NOTE -- Xhat may come out with non-zero imaginary components; if
%       your original data matrix X is strictly real valued, then you need
%       to use real(Xhat)
% 
%
% BWB, Apr 2014

%% input parsing
p = inputParser; 

% required inputs
p.addRequired('Phi', @isnumeric);
p.addRequired('lambda', @isnumeric);
p.addRequired('x0', @isnumeric);
p.addRequired('M', @(x)isnumeric(x) && x>0);

% parameter value iputs
p.addParamValue('keep_modes', [], @isnumeric);

% now parse the inputs
p.parse(Phi, lambda, x0, M, varargin{:});
inputs = p.Results;

%% if we're only using a subset of the modes in the reconstruction
if ~isempty(inputs.keep_modes),
    Phi = Phi(:, inputs.keep_modes);
    lambda = lambda(inputs.keep_modes);
end;

%% reconstruct
z0 = Phi\x0; % initial conditions
Z = zeros(length(z0), M);
for tt = 1:M,
    Z(:, tt) = z0 .* lambda.^tt;
end;

Xhat = Phi * Z;