close all;
clear;

% Parameters
Fs = 1000; % Sampling frequency
T = 1/Fs; % Sampling period
symbol_rate = 10; % Symbol rate in symbols per second
T_symbol = 1/symbol_rate; % Symbol period

% Generate test signal
t = 0:T:T*999; % Time vector
symbols = sign(randn(1, length(t))); % Random symbols (-1 or 1)
signal = symbols;%repmat(symbols, 1, 10); % Repeat symbols to create signal

% Apply different sampling phase offsets
phase_offsets = [0, 0.25, 0.5]; % Fraction of symbol period
for offset = phase_offsets
    % Create received signal with sampling phase offset
    received_signal = signal .* cos(2*pi*symbol_rate*t + 2*pi*offset);
    
    % Perform timing recovery
    [recovered_signal, timing_error] = timing_recovery_gardner(received_signal, T_symbol, T);%T, round(T_symbol/2));

    % Plot original and received signals
    figure;
    ax1 = subplot(2, 1, 1);
    stem(t, signal);
    hold on;
    stem(t, received_signal, 'LineStyle','none');
    title(['Original Signal vs. Received Signal (Phase Offset = ', num2str(offset), ')']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Original Signal', 'Received Signal');

    % Plot recovered signal
    ax2 = subplot(2, 1, 2);
    stem(t, signal);
    hold on;
    stem(t, recovered_signal, 'LineStyle','none');
    title('Original Signal vs. Recovered Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Original Signal', 'Recovered Signal');

    linkaxes([ax1, ax2], 'x');
    xlim([0 0.25]);
end
