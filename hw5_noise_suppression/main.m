clc;
clear all;
% Visually compare the results of denoising in a sig signal (file attached) using:
%   FIR filtering with a low-pass filter of the order of e.g. 40 (fir1 and filter in Matlab), 
%       estimate the cutoff frequency from the frequency content of the signal obtained by e.g. 
%       the Welch PSD estimation method (function pwelch in Matlab)
%   Filtering with a median filter (e.g. length 7, medfilt1 in Matlab)
%   DWT using the waveletSignalDenoiser tool, default options Denoising method: Bayes 
%       (empirical Bayes thresholding), Thresholding Rule: Median, Noise Estimate: Level-Independent,
%       use wavelets: db1 (Daubechies, db1=Haar wavelet) and sym2 (Symlet) and level 
%       (number of decomposition levels) 2 and 5 for both wavelets (i.e. 4 combinations in total)

load sig
% loaded data: b f fs out outhi pxx sig

% Welch PSD estimation method shows low freq of 0.1
figure(1)
subplot(2, 1, 1)
hold on
title("sig.mat pxx")
plot(pxx)

subplot(2, 1, 2)
px = pwelch(sig)
hold on
title("px = pwelch(sig)")
plot(px)

figure(2)
subplot(2, 2, 1)
hold on;
title("sig.mat out")
plot(out)

subplot(2, 2, 2)
hold on;
title("sig.mat outhi")
plot(outhi)

% FIR filtering with a low-pass filter of the order of e.g. 40
subplot(2, 2, 3)
plot(sig)
hold on;
ftr = fir1(40, 0.1, 'low');
% plot(ftr)
% freqz(ftr,1);
ot = filter(ftr, 1, sig);
xlabel("Samples")
ylabel("Amplitude")
title("low pass filter, window 40")
plot(ot)

% Filtering with a median filter
subplot(2, 2, 4)
plot(sig)
hold on;
ot = medfilt1(sig, 7);
xlabel("Samples")
ylabel("Amplitude")
title("median filter, length 8")
plot(ot)

% Wavelet toolbox
figure(3)
subplot(2, 2, 1)
sig1 = wdenoise(sig,2, ...
    Wavelet='db1', ...
    DenoisingMethod='Bayes', ...
    ThresholdRule='Median', ...
    NoiseEstimate='LevelIndependent');
hold on;
xlabel("Samples")
ylabel("Amplitude")
title("db1 2 Bayes Median LevelIndependent")
plot(sig1)

subplot(2, 2, 2)
sig2 = wdenoise(sig,5, ...
    Wavelet='db1', ...
    DenoisingMethod='Bayes', ...
    ThresholdRule='Median', ...
    NoiseEstimate='LevelIndependent');
hold on;
xlabel("Samples")
ylabel("Amplitude")
title("db1 5 Bayes Median LevelIndependent")
plot(sig2);

subplot(2, 2, 3)
sig3 = wdenoise(sig,2, ...
    Wavelet='sym2', ...
    DenoisingMethod='Bayes', ...
    ThresholdRule='Median', ...
    NoiseEstimate='LevelIndependent');
hold on;
xlabel("Samples")
ylabel("Amplitude")
title("sym2 2 Bayes Median LevelIndependent")
plot(sig3);

subplot(2, 2, 4)
sig4 = wdenoise(sig,5, ...
    Wavelet='sym2', ...
    DenoisingMethod='Bayes', ...
    ThresholdRule='Median', ...
    NoiseEstimate='LevelIndependent');
hold on;
xlabel("Samples")
ylabel("Amplitude")
title("sym2 5 Bayes Median LevelIndependent")
plot(sig4);
