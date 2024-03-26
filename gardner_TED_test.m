close all;
clear;

% Values that we're giving you
Fs = 4e6; % Sampling rate (samples/sec)
Rs = 50e3; % Symbol rate in symbols/sec (baud)
Tmax = 0.1; % Max time for the simulation (sec)
fc = 1e6; % Carrier frequency (Hz)

% Complete these expressions using the variables above
N = Fs*Tmax; % Total number of sample points in the simulation
Ns = Rs*Tmax; % Number of symbols to send
sps = N/Ns; % Number of samples per symbol.

%% 

loadStruct = load('symb_sync_input.mat');
Ybb = loadStruct.Ybb;

loadStruct = load('module4_symbols.mat');
transSymbols = loadStruct.module4_symbols;

inph_Ybb = real(Ybb);
quad_Ybb = imag(Ybb);

[recovered_signal, timing_error] = timing_recovery_gardner(Ybb, 1/Rs, 1/Fs);

%%
figure;
plot(real(Ybb));
title('Symbol sync input (after Costas loop)');
xlim([0 10000]);

scatterplot(Ybb(1:sps:N));
title('IDEAL Received constellation (sampled perfectly)');

scatterplot(recovered_signal(1:sps:N));
title('Gardner TED function');

figure;
plot(timing_error);
title('Gardner TED timing error');