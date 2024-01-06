function [tx_signal, rx_apriori] = run_tx_sim(constants)
%% Setup

% Use these variables for plotting
%% Symbol Generation
% Generate random set of bits of 0s and 1s

num_data_bits = length(constants.data_bits);
% Map bits to BPSK symbols of -1s and 1s
symbs = 2 * (constants.data_bits - 0.5);
%% Packet creation

% packet header = 3 keys

packet_header = repmat(constants.key,1,3);

packet = cat(2, packet_header, symbs);
delay = round(normrnd(100, 5)); % in samples, (=0.25ms)
num_packets = fix((constants.Ns-delay)/length(packet));
bits_remaining = rem((constants.Ns-delay),length(packet));

transmit_symbs = [zeros(1, delay) repmat(packet,1,num_packets) zeros(1,bits_remaining)];
transmit_samps = upsample(transmit_symbs, constants.sps);

figure
% stem(t,transmit_symbs, "DisplayName", "Delta train bits");
% xlim([0.001 0.003]);
% legend('show');
% Plot BPSK constellation digram
% figure
% scatterplot(symbs);
% title('BPSK Constellation Diagram');

%% Carrier Modulation

% delta_train_symbs = upsample(symbs, sps);
% figure
% plot(t, delta_train_symbs, "DisplayName", "Delta train bits");
% xlim([0 0.001]);

% rect_impulse = ones([1, constants.sps]);
roll_off =0.67;
span =10; % the number of symbols the impulse response of the RRC impulse response will last for
rrc_impulse = rcosdesign(roll_off, span, constants.sps, "sqrt");

symb_wave = conv(transmit_samps, rrc_impulse, "same");

hold on
plot(constants.t, symb_wave, "DisplayName", "RRC-filtered signal");

% Modulation
carrier = exp(-1i*2*pi*constants.fc*constants.t); % needed the 2*pi here
analytic_sig = carrier .* symb_wave;

ph_offset = 0;%unifrnd(0, 2*pi);
% select_ph_offset =0;
disp(["Phase: ", ph_offset]);
Tx = real(analytic_sig * exp(-1i*ph_offset));
hold on
plot(constants.t, Tx, 'DisplayName', 'Tx signal');
xlim([0.001 0.003]);

title('Carrier Modulation');
xlabel('Time (s)');
ylabel('Amplitude (V)');
legend('show');

tx_signal = Tx;
rx_apriori = [packet_header, rrc_impulse, transmit_samps];
end