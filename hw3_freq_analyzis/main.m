clc;
clear all;
% Create a sine signal with a frequency of 100Hz (sampling frequency 1kHz) 
% with amplitude modulation (e.g. using function 'square') so that the modulation 
% signal has a frequency of 5 Hz, a square wave with an alternation of 5% 
% and varies the amplitude of the carrier sine wave to 0.5 and 1, 
% see Fig. Do not forget to recalculate the correct scale of the axes of the resulting graphs.

% Add 0.01*randn noise to the signal

% Calculate and show:
%   signal envelope
%   immediate phase
%   instantaneous frequency
%   two-sided amplitude spectrum (using DFT)
%   two-sided spectrum of the analytical signal
%   envelope amplitude spectrum

f_sampling = 1000; % [Hz]
t_beg = -pi/64;
t_end = pi/12;       % [s]

t = t_beg : 1/f_sampling : t_end - 1/f_sampling;

fs = 100; % [Hz]
x = sin(2*pi*fs*t); % sinus with frequency 100Hz

fs_mod = 5;
A = 0.25*square(2*pi*fs_mod*t, 5) + 0.75; % amplifier
y = A.*x;

% add noise 
y = y + 0.01*randn(size(y));

% Signal itself
% figure(1);
% subplot(3, 2, 1);
% % plot(t, x);
% title("Original signal");
% hold on;
% plot(t, y);
% hold off;
% xlabel('Time [s]');
% ylabel('Values [-]');
% grid on;

% Signal envelope
subplot(3, 2, 1)
title("Original signal with envelope");
y_hat = hilbert(y);
env = abs(y_hat);
hold on;
plot(t, y);
plot(t, [-1; 1] * env, 'linewidth',1);
xlabel('Time [s]');
ylabel('Values [-]');
legend('q','up','lo');
grid on;

% Instant phase 
subplot(3, 2, 2)
phase = unwrap(angle(y_hat));
plot(t, phase)
xlabel('Time [s]');
ylabel('Phase');
grid on;
title("Instant phase");

% Instant frequency
subplot(3, 2, 3)
% instfreq(y,f_sampling, "Method", "hilbert")
instfrq = f_sampling/(2*pi)*diff(unwrap(angle(y_hat)));
plot(t(2:end),instfrq);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
grid on;
title("Instant frequency");

% Two-sided amplitude spectra
subplot(3, 2, 4)
Y = fft(y)./length(y);

Nf = length(Y);
% f = 0 : 1 : Nf - 1;
f = -Nf/2 : 1 : Nf/2 - 1;
f = f_sampling/Nf*f;

plot(f, fftshift(abs(Y)));
hold on
xlim([-f_sampling/2 f_sampling/2])
xlabel('Original signal spectrum [Hz]');
ylabel('|Y| [-]');
title("Two-sided amplitude spectra");

% =====================================================
% Two-sided spectra of analytic signal
subplot(3, 2, 5)
Y = fft(y_hat)./length(y);

Nf = length(Y);
f = -Nf/2 : 1 : Nf/2 - 1;
f = f_sampling/Nf*f;

plot(f, fftshift(Y));
xlim([-f_sampling/2 f_sampling/2])
hold on;
xlabel('Analytic signal spectrum [Hz]');
% xlim([-f_sampling f_sampling])
ylabel('|Y| [-]');
title("Two-sided spectra of analytic signal");

% Amplitude envelope spectra
subplot(3, 2, 6)
Y = fft(real(env))./length(y);

Nf = length(Y);
f = 0 : 1 : Nf - 1;
f = f_sampling/Nf*f;


plot(f, fftshift(abs(Y)));
xlabel('Envelope spectrum [Hz]');
ylabel('|Y| [-]');
title("Amplitude envelope spectra");

