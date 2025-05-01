%% Simulation Setup

Fc = 1e6; % Carrier frequency in Hz
Fs = 4e6; % Sampling rate
Rs = 50e3; % Symbol rate in symbols/sec (baud)
B = 350e3; % Bandwidth

sps = Fs / Rs; % Number of samples per symbol
num_symbols = 368; % Number of symbols/bits to send
num_samples = num_symbols * sps; % Number of discrete samples

% Use these guys for plotting
t = linspace(0, num_symbols / Rs, num_samples); % Time variable
f = linspace(-Fs/2,Fs/2, num_samples); % Frequency variable

% Clear all previous figures before starting
close all;

%% Mapping Bits to Symbols

% Below are the bits for the binary message that we want to send

bits = [0 1 0 0 1 0 0 0 0 1 1 0 0 0 0 1 0 1 1 1 0 0 1 0 0 1 1 0 0 1 0 0 0 1,... 
        1 1 0 1 1 1 0 1 1 0 0 0 0 1 0 1 1 1 0 0 1 0 0 1 1 0 0 1 0 1 0 0 1 0,...
        0 0 0 0 0 1 0 1 0 0 1 1 0 1 1 1 0 1 0 1 0 1 1 0 0 0 1 1 0 1 1 0 1 0,...
        1 1 0 1 1 1 0 0 1 1 0 0 1 0 0 0 0 0 0 1 0 1 0 0 1 1 0 1 1 1 0 1 1 1,...
        0 1 1 0 1 0 0 1 0 1 1 1 0 1 0 0 0 1 1 0 0 0 1 1 0 1 1 0 1 0 0 0 0 0,...
        1 0 0 0 0 0 0 1 0 1 0 1 0 0 0 1 1 0 1 1 1 1 0 0 1 0 0 0 0 0 0 1 0 0,...
        0 0 1 1 0 1 1 0 1 1 1 1 0 1 1 0 1 1 0 1 0 1 1 1 0 0 0 0 0 1 1 1 0 1,...
        0 1 0 1 1 1 0 1 0 0 0 1 1 0 0 1 0 1 0 1 1 1 0 0 1 0 0 0 1 0 0 0 0 0,...
        0 1 0 0 0 1 0 1 0 1 1 0 1 1 1 0 0 1 1 0 0 1 1 1 0 1 1 0 1 0 0 1 0 1,...
        1 0 1 1 1 0 0 1 1 0 0 1 0 1 0 1 1 0 0 1 0 1 0 1 1 1 0 0 1 0 0 1 1 0,...
        1 0 0 1 0 1 1 0 1 1 1 0 0 1 1 0 0 1 1 1 0 0 1 0 0 0 0 1];

% Convert to BPSK Symbols (1's and -1's)
symbols = bits;
symbols(symbols == 0) = -1;

% Visualize BPSK Symbols in a Constellation Diagram
scatterplot(symbols);
title("Transmitted Symbols");

%% Represent Symbols as Time Shifted Deltas

% Use upsample to form a delta train
deltas = upsample(symbols, sps);

% Visualize Deltas (xlimited to only first ~20 symbols)
figure;
stem(t, deltas);
title("Upsampled Deltas");
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 0.0004]);

%% Filter Deltas with Transmitter Filter

% Use an SRRC as the pulse shape filter 
beta = 0.5; % RRC rolloff factor
span = 5; % number of symbols for length of filter impulse response
ps_filter = rcosdesign(beta, span, sps, 'sqrt'); 

% Convolve the deltas with the rectangular window
transmited_baseband = conv(deltas, ps_filter, 'same'); 

% Visualize Transmitted Baseband Signal
figure;
plot(t, transmited_baseband);
title("Transmited Baseband Signal");
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 0.0004]);

%% Modulate Baseband Signal to Passband

% TODO 1.1: Edit modulate_carrier to implement a constant phase_offset
% Hint: Add phase_offset as an input paramter to modulate_carrier

phase_offset = pi/4;
transmitted_signal = modulate_carrier(transmited_baseband, Fc, t, phase_offset);

%% Transmit Through Wireless Channel

snr = 10;
received_signal = add_channel_impairments(transmitted_signal, B, Fc, Fs, snr);

% Plot the Received Signal
figure;
plot(t, received_signal)
title("Received Signal")
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 0.0004]);

%% Demodulation

% TODO 1.2: Encorporate a frequency offset
frequency_offset = Fc * 0.05;
Fr = Fc + frequency_offset;

% TODO 3.0.0 Comment below line out so we can replace with the costas_loop function
% [I, Q] = naive_demod(received_signal, Fr, Fs, t);

% TODO 3.0.1: and normalize the signal for stability of PID tuning constants
% Hint: Use the normalize function
normalized_signal = normalize(received_signal);

% PID Tuning Constants
Kp = 0.1;
Ki = 0.002;
Kd = 0;

% TODO: Uncomment line below and complete Costas Loop function
[I, Q, theta, err] = costas_loop(normalized_signal, Fr, Fs, t, Kp, Ki, Kd);

figure;
plot(t, I)
title("Received In Phase Component")
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 0.0004]);

%% Apply Receiver Impulse Response

% Combine I and Q components to construct the received baseband signal (I + jQ)
received_baseband = I + 1j * Q;

% Convolve the received baseband signal with the receiver filter
received_samples = conv(received_baseband, ps_filter, 'same'); 

%% Sample and Detect Symbols

% (Naively) sample the received signal at the symbol rate to get the received symbols
% Fixing this assumption will be a main focus of the next module
received_symbols = downsample(received_samples, sps);

% Visualize the received symbols in a constellation diagram (scatterplot)
scatterplot(received_symbols);
title("Received Symbols");

% Apply Thresholding to detect the symbols (turning them back to 1's and -1's)
received_symbols(received_symbols > 0) = 1;
received_symbols(received_symbols <= 0) = -1;

%% Map Symbols back to Bits

% Map the detected symbols back to bits (turn them back to 1's and 0's
detected_bits = received_symbols;
detected_bits(detected_bits == -1) = 0;

% Calculate Bit Error Rate
num_bit_errors = sum(detected_bits(1:length(bits)) ~= bits)
BER = 100 * num_bit_errors / length(bits)

%% Decode and Display Message

% If you have successfully decoded the message, you should see a nice
% message
s = string(detected_bits);
s = strjoin(s);
s = strrep(s, " ", "");
disp(s);
s_len = length(char(s))
inputString = char(s);
binaryString = inputString(1:end-mod(length(inputString),8));
binaryChunks = reshape(binaryString, 8, []).';
asciiChars = char(bin2dec(binaryChunks)).';
disp(asciiChars);
