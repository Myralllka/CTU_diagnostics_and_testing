clc;
clear all;
syms x y z xi yi zi xj yj zj v dij;

eq = sqrt((x-xj)^2 + (y-yj)^2 + (z-zj)^2) - sqrt((x-xi)^2 + (y-yi)^2 + (z-zi)^2) - v*dij;
% diff(eq, x)
simplify(diff(eq, y))
% diff(eq, y)
% diff(eq, z)