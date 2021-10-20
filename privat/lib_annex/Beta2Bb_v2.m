% Beta2Bb.m
% kniewiad@obs-vlfr.fr  2006-07-17  Initial
function [beta_p, Bbp] = Beta2Bb_v2(S, lambda, beta)
% This script calculates particulate backscatter at a chosen wavelength from the
% scattering

% Inputs
% S is mean salinity of data set
% lambda is scatter wavelength of sensor
% Beta is the scatter at lambda

X = 1.1; % See Boss and Pegau, 2001 and Wetlabs manual (1.077 & 1.132 for Eco-BB & FLNTU)
% Wetlabs ECOBB3 Manual
delta = 0.09;
theta = 124; % Scatter angle (124 for Wetlabs Eco-BB and 140 for FLNTU from Sulivan et al., 2010)
% Assuming a linear relationship with salinity and using table 4 of Morel 1974
A = 1.38*(lambda/500)^-4.32*(1+0.3*S/37)*10^-4;

beta_w = A*(((1+cos(theta)^2)*(1-delta))/(1+delta));  % Morel 1974;

beta_p = beta - beta_w;

% Calculated backscatter from scatter of ECO puck
Bbp=2*pi*beta_p*X;
end