function [Phi, mu, lambda, diagS, x0] = DMD(Xraw, varargin)
% function [Phi mu lambda diagS x0] = DMD(Xraw, varargin)
% computes the Dynamic Mode Decomposition of data matrix Xraw
%
% INPUTS: 
%
% rows of Xraw are assumed to be measurements
% columns Xraw are assumed to be time points, sampled at equal dt's
%
% optional parameters: {'parameter_name', [default_value]}
%   {'dt', [1]}
%   {'r', [1e32]}        truncate svd basis to first r features
%   {'scale_modes', 1}   0 or 1, scale dmd modes by singular values (energy)
%   {'nstacks', 1}       number of stacks of the raw data

%
% OUTPUTS:
%
% Phi, the modes
% mu, the fourier spectrum of modes (mu = log(lambda)/dt)
% lambda, the DMD spectrum of modes 
% diagS, the singular values of the data matrix
% x0, the initial condition vector corresponding to Phi
% 
%
% BWB, Nov 2013
% mods: return diag(S), BWB Dec 2013
%       added corrected truncation of modes, BWB Jan 2014
%       added scaled modes, BWB Mar 2014
%       added stacking to get augmented data matrix, BWB Mar 2014
%       added option to use optimal singular value hard threshold (SVHT,
%             see Gavish & Donoho 2013), BWB Jun 2014

%% input parsing
p = inputParser; 

% required inputs
p.addRequired('Xraw', @isnumeric);

% parameter value iputs
p.addParameter('dt', 1, @(x)isnumeric(x) && x>0);
p.addParameter('r', 1e32, @(x)isnumeric(x) && x>0);
p.addParameter('use_optimal_SVHT', 0, @isnumeric);
p.addParameter('scale_modes', 1, @isnumeric);
p.addParameter('nstacks', 1, @(x)isnumeric(x) && x>0);

% now parse the inputs
p.parse(Xraw, varargin{:});
inputs = p.Results;

%% stacking the data matrix 
if inputs.nstacks > 1,
    Xaug = [];
    for st = 1:inputs.nstacks,
        Xaug = [Xaug; Xraw(:, st:end-inputs.nstacks+st)]; %#ok<AGROW>
    end;
    
    X = Xaug(:, 1:end-1);
    Y = Xaug(:, 2:end);
else
    X = Xraw(:, 1:end-1);
    Y = Xraw(:, 2:end);
end;

%% DMD
[U, S, V] = svd(X, 'econ');
diagS = diag(S);

% if we want to use optimal singular value hard threshold, compute the
% truncation order r
if inputs.use_optimal_SVHT > 0,
    beta = size(X,1)/size(X,2); if beta > 1, beta = 1/beta; end;
    omega = optimal_SVHT_coef(beta,0) * median(diagS);
    r = sum(diagS > omega);
else
    r = inputs.r;
end;

if r >= size(U,2),
    % no truncation
    Atilde = U'*Y*V/S;
    
    if inputs.scale_modes == 0,
        [W, D] = eig(Atilde);  
    else % scaling modes
        Ahat = S^(-1/2) * Atilde * S^(1/2);
        [What, D] = eig(Ahat);
        W = S^(1/2) * What;
    end;

    Phi = Y*V/S*W;
else
    % truncate modes
    U_r = U(:, 1:r);
    S_r = S(1:r, 1:r);
    V_r = V(:, 1:r);
    
    Atilde = U_r' * Y * V_r / S_r;
    
    if inputs.scale_modes == 0,
        [W_r, D] = eig(Atilde);
    else % scaling modes
        Ahat = (S_r^(-1/2)) * Atilde * (S_r^(1/2));
        [What, D] = eig(Ahat);
        W_r = S_r^(1/2) * What;
    end;
    
    Phi = Y*V_r/S_r*W_r;
end;

lambda = diag(D);
mu = log(lambda)/inputs.dt;
x0 = X(:,1);
 