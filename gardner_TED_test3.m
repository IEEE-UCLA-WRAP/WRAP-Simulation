close all;
clear;

% Parameters
fs = 4e6; % Sampling frequency
T = 1/fs; % Sampling period
samples_per_symbol = 80; % Upsampling factor
Rs = T*samples_per_symbol;
symbols = randi([0 1], 1, 1000); % Random binary data
modulated_symbols = 2*symbols - 1; % BPSK modulation

% Upsample at the transmitter
tx_signal = upsample(modulated_symbols, samples_per_symbol);

% Add noise
SNR_dB = 10; % Signal-to-Noise Ratio in dB
noise_power = 10^(-SNR_dB/10);
noise = sqrt(noise_power) * randn(1, length(tx_signal));
rx_signal = tx_signal + noise;

% Apply different sampling phase offsets
phase_offsets = [0, 0.25, 0.5]; % Fraction of symbol period
for offset = phase_offsets
    % Create received signal with sampling phase offset
    received_signal = tx_signal .* cos(2*pi*Rs*t + 2*pi*offset);
    
    % Perform timing recovery
    [recovered_signal, timing_error] = timing_recovery_gardner(rx_signal, Rs, T);%T, round(T_symbol/2));

    % Plot original and received signals
    figure;
    ax1 = subplot(2, 1, 1);
    stem(t, tx_signal);
    hold on;
    plot(t, rx_signal);%stem(t, received_signal, 'LineStyle','none');
    title(['Original Signal vs. Received Signal (Phase Offset = ', num2str(offset), ')']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Original Signal', 'Received Signal');

    % Plot recovered signal
    ax2 = subplot(2, 1, 2);
    stem(t, signal);
    hold on;
    plot(t, recovered_signal);%stem(t, recovered_signal, 'LineStyle','none');
    title('Original Signal vs. Recovered Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Original Signal', 'Recovered Signal');

    linkaxes([ax1, ax2], 'x');
    xlim([0 0.001]);

    % Perform demodulation
    demodulated_symbols = recovered_signal(1:samples_per_symbol:end) > 0;
    
    % Calculate bit error rate (BER)
    ber = sum(symbols ~= demodulated_symbols) / length(symbols);

    fprintf('Bit Error Rate (BER): %.5f\n', ber);

end

% % Plot results
% figure;
% subplot(2,1,1);
% plot(rx_signal);
% title('Received Signal with Noise and Timing Offset');
% xlabel('Sample');
% ylabel('Amplitude');
% 
% subplot(2,1,2);
% stem(timing_error);
% title('Timing Error');
% xlabel('Sample');
% ylabel('Timing Error (T)');
