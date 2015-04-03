function [f P] = DMD_spectrum(Phi, mu, varargin)
% function [f P] = DMD_spectrum(Phi, mu, varargin)
% computes the DMD spectrum, based on outputs from DMD.m
%
% INPUTS: 
%
% Phi and mu are outputs from DMD.m
%
% optional parameters: {'parameter_name', [default_value]}
%   {'plotit', [0]}        0 or 1, whether to plot the spectrum
%
% OUTPUTS:
%
% f, the frequencies of the modes in cycles/sec
% P, the power of the modes
%
% BWB, Apr 2014

%% input parsing
p = inputParser; 

% required inputs
p.addRequired('Phi', @isnumeric);
p.addRequired('mu', @isnumeric);

% parameter value iputs
p.addParamValue('plotit', 0, @isnumeric);

% now parse the inputs
p.parse(Phi, mu, varargin{:});
inputs = p.Results;

%%
f = abs(imag(mu(:))/2/pi); % frequency in cycles/sec
P = (diag(Phi'*Phi)); % roughly scales like the fft spectrum

f = f(:);
P = P(:);
%%
if inputs.plotit,
    figure;
    stem(f, P, 'k');
    axis tight;
    xlabel('freq. (Hz)');
    ylabel('mode power');
end;