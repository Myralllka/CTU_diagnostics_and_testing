clc;
clear all;
syms x y z xi yi zi xj yj zj v dij;

% eq = abs(sqrt((x-xj)^2 + (y-zyj)^2 + (z-zj)^2) - sqrt((x-xi)^2 + (y-yi)^2 + (z-zi)^2)) - v*dij;
% diff(eq, x)
% simplify(diff(eq, x))
% diff(eq, y)
% diff(eq, z)
% =====================================
 dw = 1/(499.2*128 * 10^6); % [s]
 c = 299792458; % [m/s]
 c_dw = c * dw; % [m/dw]
% ==================================

coors = [-1.97, -12.75, -12.77, -1.81, -6.86, -1.92, -6.87, -12.27, -6.77;
        -8.05, -8.05, 2.75, 2.75, -2.67, -2.67, -8.05, -2.67, 2.75;
        2.6, 2.6, 3.13, 3.13, 2.86, 2.86, 2.6, 2.86, 3.13];

corrections = [25882948.1690979,
               204304317.49362183,
               -852505246.5370178,
               -312401912.46502686,
               0,
               -113743765.83856201,
               87721032.96951294,
               -25450735.08444214,
               -1209247150.8052368];
TOF = zeros(9, 9);
for i = 1:9
    for j = 1:9
        TOF(i, j) = norm(coors(:, i) - coors(:, j)) / c_dw;
    end
end
disp(TOF)

M = csvread("syncs_cold_start.csv");
% M(:, 5) = M(:, 5) - M(1, 5);
% a(a(:, 1) == 1, :)
M = M(8200:end, :);
s1 = M(M(:, 1) == 68, :);
s2 = M(M(:, 1) == 69, :);
s3 = M(M(:, 1) == 70, :);
s4 = M(M(:, 1) == 71, :);
s5 = M(M(:, 1) == 72, :);
s6 = M(M(:, 1) == 73, :);

% assume s1 as main.

a = 1;
b = 1000;

hold on
% plot(s1(a:b, 5), s1(a:b, 3) + corrections(1));
% plot(s2(a:b, 5), s2(a:b, 3) + corrections(2));
% plot(s3(a:b, 5), s3(a:b, 3) + corrections(3));
plot(s4(a:b, 5), s4(a:b, 3) + corrections(4) - TOF(5, 4));
plot(s5(a:b, 5), s5(a:b, 3) + corrections(5));
% scatter(s6(a:b, 5), s6(a:b, 3));
