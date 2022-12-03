clc;
clear all;
syms x y z xi yi zi xj yj zj v dij;

eq = abs(sqrt((x-xj)^2 + (y-zyj)^2 + (z-zj)^2) - sqrt((x-xi)^2 + (y-yi)^2 + (z-zi)^2)) - v*dij;
% diff(eq, x)
simplify(diff(eq, x))
% diff(eq, y)
% diff(eq, z)

% dw = 1/(499.2*128 * 10^6); % [s]

% c = 299792458; % [m/s]

% speed of light in m/dw:
% cdw = c / dw