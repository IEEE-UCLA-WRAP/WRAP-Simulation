close all

% Values that we're giving you
Fs = 5e6; % Sampling rate (samples/sec)
Rs = 40e3; % Symbol rate in symbols/sec (baud)
Tmax = 0.1; % Max time for the simulation (sec)
fc = 1e6; % Carrier frequency (Hz)

% Complete these expressions using the variables above
N = 25600; % Total number of sample points in the simulation
Ns = 256; % Number of symbols to send
sps = 100; % Number of samples per symbol.

% Use these variables for plotting
t = linspace(0, Ns / Rs, N + 1); % Time vector

% Generate the analytic signal by multiplying the convolved signal with exp(-j * 2 * pi * fc * t)
cc = cos( 2 * pi * fc * t);

figure;
plot(real(cc));

min(real(cc))

