%% OFDM_sandbox.m
%
% Author: Abraham Canafe
%
% Adapted from https://dspillustrations.com/pages/posts/misc/python-ofdm-example.html
%
% Date: 4/17/2024

clearvars, close all, clc
rng(42)
%% Parameters
K = 1024; % number of OFDM subcarriers (ensure that this is a power of 2)
CP = floor(K/4); % cyclic prefix size
QAM_mod_idx_M = 16;

% Define parameters for pilot symbols
P = 64; % number of pilot symbols; increasing pilot symbols results in better channel estimates
pilot_value_integer = 8; % (1100 in bits)
pilot_value = qammod(pilot_value_integer, QAM_mod_idx_M,UnitAveragePower=true); % QAM representation

all_carriers = (1:K)';

pilot_carriers = all_carriers(1:floor(K/P):K);
pilot_carriers = [pilot_carriers; all_carriers(end)];

data_carrier_bools = ~ismember(all_carriers,pilot_carriers);
data_carriers = all_carriers(data_carrier_bools);

figure(1)
plot(data_carriers,zeros(size(data_carriers)),'ob')
hold on
plot(pilot_carriers,zeros(size(pilot_carriers)),'or')
legend('Data subcarriers','Pilot subcarriers','interpreter','latex')
xlabel('Carrier index','interpreter','latex')
hold off


%% Channel response
% ============== Play around with the channel size/type ============== %
% We normalize channel gains for convenience.
% Note that the magnitude of a complex Gaussian follows a Rayleigh
% distribution and that the angle of a single-tap Gaussian channel follows
% a uniform distribution.

% ===== Ideal channel ===== %
% h = 1;
% go figure

% ===== Rayleigh fading ===== %
%h = 1/sqrt(2) * (randn(1) + 1j.*randn(1));

% ===== Arbitrary channel from Python tutorial ==== %
h = [1, 0, 0, 0, 0, 0.3 + 0.3j] ./ norm([1, 0, 0, 0, 0, 0.3 + 0.3j]);
% try increasing/decreasing the number of pilots to estimate this channel

% ===== Multiple Gaussian taps ===== %
% n_taps = 1; % number of delay taps for channel size > 1
% channel_size = [1,n_taps]; 
% h = 1/sqrt(2) * (randn(channel_taps) + 1j.*randn(channel_taps));

% =================================================================== %
H_fft = fft(h,K);
figure(2)
plot(all_carriers,abs(H_fft));
xlabel("Subcarrier Index",'interpreter','latex')
ylabel("$|H(f)|$",'interpreter','latex');
title("Channel impulse response",'interpreter','latex')

%% Generate data
data_bitstream = randi([0 1], length(data_carriers)*log2(QAM_mod_idx_M), 1);

%% Symbol mapping
data_qam_syms = qammod(data_bitstream, QAM_mod_idx_M, InputType='bit', UnitAveragePower=true);
pilot_qam_seq = repmat(pilot_value, length(pilot_carriers), 1);
scatterplot(data_qam_syms)
title('Tx constellation')

%% Map QAM symbols to appropriate subcarriers
OFDM_data = zeros(K,1);
OFDM_data(pilot_carriers) = pilot_qam_seq;
OFDM_data(data_carriers) = data_qam_syms;

% ========== Experimental ========== %
% Technically, when symbols are being packed into the OFDM subcarriers
% some type of "pulse shaping" in the frequency domain is supposed to occur.
% This isn't usually necessary for OFDM simulations. Of course, one can
% design the sinc pulses to ensure nulls at desired subcarrier indices.

% freq_bins = -length(OFDM_data)/2:length(OFDM_data)/2;
% freq_pulse_ker = sinc(freq_bins);
% OFDM_data = conv(OFDM_data,freq_pulse_ker,'same');
% ================================== %

%% IFFT
OFDM_data_ifft = ifft(OFDM_data);

%% Add cyclic Prefix
cp_OFDM_data_ifft = [OFDM_data_ifft(end-CP+1:end); OFDM_data_ifft];

%% Channel convolution + noise
% We assume baseband model, so we can ignore the carrier frequency/sampling rate for now.
% To include carrier frequency, sampling rate, etc. we include the following:
% conv_cp_OFDM_ifft = conv(cp_OFDM_data_ifft, h) * exp(1j*2*pi*delta_f*(0:length(cp_OFDM_data_ifft))*Ts)
% delta_f = ___; % carrier frequency
% fs = ___; % sampling rate
% Ts = 1 / fs;
conv_cp_OFDM_ifft = conv(cp_OFDM_data_ifft, h);

% ======== Add Additive White Gaussian Noise ========= %
P_cp_OFDM_ifft = mean(abs(conv_cp_OFDM_ifft.^2)); % compute signal power
SNR_dB = 20; % toggle signal-to-noise ratio
SNR_linear = 10^(SNR_dB / 10);
sigma_awgn = P_cp_OFDM_ifft / SNR_linear; % compute noise power
awgn_noise = sqrt(sigma_awgn/2) .* randn(size(conv_cp_OFDM_ifft));
received_cp_OFDM_ifft = conv_cp_OFDM_ifft + awgn_noise;
% The above is equivalent to using the single-line MATLAB function:
% received_cp_OFDM_ifft = awgn(conv_cp_OFDM_ifft, SNR_dB, 'measured');
% ==================================================== %


%% Remove CP
received_OFDM_ifft = received_cp_OFDM_ifft(CP+1:CP+K); % ifft
figure
plot(abs(received_OFDM_ifft))
hold on
plot(abs(OFDM_data_ifft))
title('Received signal with removed cyclic prefix','interpreter','latex')
xlabel('$t$','interpreter','latex')
ylabel('$|y(t)|$','interpreter','latex')
legend('Tx signal','Rx signal','interpreter','latex')
hold off

received_OFDM_fft = fft(received_OFDM_ifft); % fft

%% Channel estimation 
H_est = estimateChannel(received_OFDM_fft, pilot_carriers, pilot_qam_seq);

% Plot of actual channel magnitude vs. estimated channel magnitude
figure(2)
hold on
plot(all_carriers,abs(H_est))
legend('Actual channel','Channel estimate','interpreter','latex')
hold off

%% Channel equalization
% TLDR: We want to "undo" the effects that the channel had upon the gain
% and phase of the Rx signal (to the extent that it is possible).
rx_equalized_OFDM_fft = received_OFDM_fft ./ H_est;
rx_equalized_qam_syms = rx_equalized_OFDM_fft(data_carriers);
scatterplot(rx_equalized_qam_syms)
title('Equalized Rx constellation')

%% Demodulation 
rx_data_bitstream = qamdemod(rx_equalized_qam_syms,QAM_mod_idx_M,OutputType='bit',UnitAveragePower=true);

%% Bit Error Rate
BER = mean(data_bitstream == rx_data_bitstream)

%% Potential Tasks
% Observe symbol constellations can change as SNR is increased/decreased
% Observe how effects of SNR and/or channel affect BER
% Etc.

%% Functions
function H_est = estimateChannel(OFDM_demod, pilot_idx, true_pilot_sequence)
    pilot_syms = OFDM_demod(pilot_idx);
    H_est_at_pilots = pilot_syms ./ true_pilot_sequence;
    all_OFDM_demod_idx = transpose(1:length(OFDM_demod));
    
    H_est_abs = interp1(pilot_idx, abs(H_est_at_pilots), all_OFDM_demod_idx, 'linear');
    H_est_angle = interp1(pilot_idx, angle(H_est_at_pilots), all_OFDM_demod_idx, 'linear');
    H_est = H_est_abs .* exp(1j.*H_est_angle); % slight spaghetti transpose
end

