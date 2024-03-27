close all;
clear;

% Parameters
Fs = 4e6; % Sampling frequency
T = 1/Fs; % Sampling period
Rs = 50e3; % Symbol rate in symbols per second
T_symbol = 1/Rs; % Symbol period
Tmax = 0.1;
N = round(Fs*Tmax); % Total number of sample points in the simulation
Ns = round(Rs*Tmax); % Number of symbols to send
sps = round(N/Ns); % Number of samples per symbol.

% Generate test signal
t = linspace(0, Tmax, N); % Time vector
% symbols = sign(randn(1, round(length(t)/sps))); % Random symbols (-1 or 1)
% signal = upsample(symbols, sps);%repmat(symbols, 1, 10); % Repeat symbols to create signal
loadStruct = load('symb_sync_input.mat');
signal = loadStruct.Ybb;

loadStruct = load('module4_symbols.mat');
symbols = loadStruct.module4_symbols;

% Apply different sampling phase offsets
phase_offsets = [0, 20, 30]; % Fraction of symbol period 
% 0.5 case can't be fixed because that's the job of frame sync?
% anything above 0.25 can't be fixed because that's the job of frame sync?
for offset = phase_offsets
    % Create received signal with sampling phase offset
    % received_signal = signal .* cos(2*pi*offset); % 2*pi*Rs*t + 
    received_signal = signal(1+offset:end);% Ybb = Ybb(1+sampPhaseOffset:end);
    
    % Perform timing recovery
    [recovered_signal, timing_error] = timing_recovery_gardner(received_signal, T_symbol, T);%T, round(T_symbol/2));

    % Plot original and received signals
    figure;
    ax1 = subplot(3, 1, 1);
    plot(t, real(signal));
    hold on;
    plot(t(1:length(received_signal)), real(received_signal));%stem(t, received_signal, 'LineStyle','none');
    title(['Original Samples vs. Received Samples (Phase Offset = ', num2str(offset), ')']);
    subtitle('Shows that the offset happened');
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Original Signal', 'Received Signal');

    % Plot recovered signal
    ax2 = subplot(3, 1, 2);
    plot(t, real(signal));
    hold on;
    plot(t(1:length(received_signal)), real(recovered_signal));%stem(t, recovered_signal, 'LineStyle','none');
    title('Original Samples vs. Recovered Samples');
    subtitle('Shows that we fixed the offset');
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Original Samples', 'Recovered Samples');

    linkaxes([ax1, ax2], 'x');
    xlim([0 0.001]);

    ax3 = subplot(3, 1, 3);
    plot(real(timing_error));
    title('Timing Error');
    xlabel('Sample');
    ylabel('Timing Error (T)');

    % Perform demodulation
    demodulated_symbols = sign(real(recovered_signal(1:sps:end)));
    
    % Calculate bit error rate (BER)
    ber = sum(symbols ~= demodulated_symbols) / length(symbols);

    fprintf('Bit Error Rate (BER): %.5f\n', ber);
end
