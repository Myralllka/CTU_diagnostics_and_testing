
Ford = n * r / 120 * (1 - Bd / Pd * cos(phi));
Fird = n * r / 120 * (1 + Bd / Pd * cos(phi));
Fbd = n * r / 120 * (1 - (Bd / Pd * cos(phi))^2);
Fc = r / 120 * (1 - Bd / Pd * cos(phi));

% Measurements duration
MEASUREMENT_PERIOD = 1;
% 
% [acc_session] = init_acc();
% 
% [data,time,rate] = get_acc(acc_session, MEASUREMENT_PERIOD); 

myfft(data(:, 1), rate);

function myfft(samples, Fs)
    % sampling period
    
    % TODO: replace this with some not random initialization of points
    % points = ...
    
    % For testing now, just generate points as a sum of 2 sinusoids

    [U, f] = fft_hertz(samples, Fs);

    plot(samples);

    figure


    plot(f,U);
    xlim([0, 200]);
    title("Single-Sided Amplitude Spectrum of X(t)")
    xlabel("f (Hz)")
    ylabel("|P1(f)|")

%     hold on;


    n = 12;
    Bd = 6.35; % [mm]
    Pd = 33.15; % [mm]
    phi = 11.3 * pi / 180; %[deg]
    
    % RPM 
    r = 2050;


    Ford = n * r / 120 * (1 - Bd / Pd * cos(phi));
    Fird = n * r / 120 * (1 + Bd / Pd * cos(phi));
    Fbd = n * r / 120 * (1 - (Bd / Pd * cos(phi))^2);
    Fc = r / 120 * (1 - Bd / Pd * cos(phi));


    xline(Ford, "-r");
    xline(Fird, '-b');
    xline(Fbd, '--g');
    xline(Fc, '-y');

    xline(r / 60, '-b');
    xline(2 * r / 60, '-b');
    xline(3 * r / 60, '-b');


end


function [U,f] = fft_hertz(u, Fs) 
    L = numel(u);
    U = abs(fft(u)/L);
    U = 2*U(1:L/2+1);
    f = Fs*(0:(L/2))/L;
end

