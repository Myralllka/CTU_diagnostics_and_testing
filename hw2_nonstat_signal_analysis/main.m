clc;
clear all;
% Create a signal - a logarithmic chirp (eg function 'chirp' in Matlab) - 
% 10 seconds long, with an initial frequency of 2Hz, an end frequency of 
% 300Hz, sampled at a frequency of 1kHz. Modulate its amplitude with a 
% Gaussian window (e.g. func gausswin):


% Calculate and display the result:
%   short-time Fourier transforms - Short Time Fourier Transform, (with default function stft parameters: Hann window length 128 samples, 75% overlap, 128 DFT samples)
%   repeat with a window of length 256 samples (and 256 DFT samples)
%   Hilbert-Huang transforms (functions emd and hht)
%   display the result of point 3 in 3D using function mesh (the result can be improved by the options 'EdgeColor','none','FaceColor','interp')

f_sampling = 1000;
t_beg = 0;
t_end = 10;
t = t_beg : 1/f_sampling : t_end - 1/f_sampling;

f_chirp_beg = 2;
f_chirp_end = 300;
x = chirp(t, f_chirp_beg, f_sampling/10, f_chirp_end);

% fs_mod = 5;
A = gausswin(length(x)); % amplifier
y = A.'.*x;

% Signal itself
figure(1);
subplot(3, 2, [1, 2])
title("Original signal");
hold on;
plot(t, y);
hold off;
xlabel('Time [s]');
ylabel('Values [-]');
title("Signal itself")
grid on;

% Short time fourier transform, default settings
subplot(3, 2, 3)
stft(y, f_sampling)
title("STFT, default settings")

% STFT, window size 256, Ndft samples 256
subplot(3, 2, 4)
M = 256;
g = hann(M,"periodic");
L = 0.75*M;
Ndft=M;
stft(x,f_sampling, Window=g,OverlapLength=L,FFTLength=256);
title("STFT, window size 256, Ndft samples 256")

% Hilbert-Huang transform
subplot(3, 2, 5);
imf = emd(y);
hht(imf,'FrequencyLimits',[0 pi]);

% Hilbert-Huang transform in 3D
subplot(3, 2, 6);
imf = emd(y, "Display", 1);
[hs,f,t] = hht(imf, 'FrequencyLimits',[0 pi]);
mesh(seconds(t),f,hs,'EdgeColor','none','FaceColor','interp')
title("aP")
xlabel('Time (s)')
ylabel('Frequency (Hz)')
zlabel('Instantaneous Energy')
