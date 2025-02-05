%% Simulation Setup

Fc = 1e6; % Carrier frequency in Hz
Fs = 4e6; % Sampling rate
Rs = 50e3; % Symbol rate in symbols/sec (baud)
B = 350e3; % Bandwidth

sps = Fs / Rs; % Number of samples per symbol
num_symbols = 460; % Number of symbols/bits to send
num_samples = num_symbols * sps; % Number of discrete samples

% Use these guys for plotting
t = linspace(0, num_symbols / Rs, num_samples); % Time variable
f = linspace(-Fs/2,Fs/2, num_samples); % Frequency variable

% Clear all previous figures before starting
close all;

%% Mapping Bits to Symbols

% Below are the bits for the binary message that we want to send

bits = [0 1 0 0 0 1 0 1 1 1 0 0 1 0 0 1 1 0 0 1 0 0 0 1 1 1 0 1 1 1 0 1 1 0,...
        0 0 0 1 0 1 1 1 0 0 1 0 0 1 1 0 0 1 0 1 0 0 1 0 0 0 0 0 0 1 0 1 0 0,...
        1 1 0 1 1 1 0 1 0 1 0 1 1 0 0 0 1 1 0 1 1 0 1 0 1 1 0 1 1 1 0 0 1 1,...
        0 0 1 0 0 0 0 0 0 1 0 1 0 0 1 1 0 1 1 1 0 1 1 1 0 1 1 0 1 0 0 1 0 1,...
        1 1 0 1 0 0 0 1 1 0 0 0 1 1 0 1 1 0 1 0 0 0 0 0 0 0 1 1 0 1 1 0 1 1,...
        1 1 0 1 1 0 1 1 0 1 0 1 1 1 0 0 0 0 0 1 1 1 0 1 0 1 0 1 1 1 0 1 0 0,...
        0 1 1 0 0 1 0 1 0 1 1 1 0 0 1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 1 0 1 0 1,...
        1 0 1 1 1 0 0 1 1 0 0 1 1 1 0 1 1 0 1 0 0 1 0 1 1 0 1 1 1 0 0 1 1 0,...
        0 1 0 1 0 1 1 0 0 1 0 1 0 1 1 1 0 0 1 0 0 1 1 0 1 0 0 1 0 1 1 0 1 1,...
        1 0 1 0 1 0 0 1 0 1 1 0 1 0 1 1 1 1 1 1 1 1 0 0 0 1 0 0 1 0 0 1 0 1,...
        0 0 1 1 0 1 1 0 1 0 1 1 0 1 1 0 1 0 0 1 0 1 1 0 0 0 1 0 0 1 1 0 1 0,... 
        0 1 0 1 1 0 0 1 0 0 0 1 1 0 1 0 0 1 0 0 1 1 0 1 1 0 1 0 1 1 0 1 1 0,... 
        1 0 0 1 0 1 1 0 0 0 1 0 0 1 1 0 1 0 0 1 0 1 1 0 0 1 0 0 0 1 1 0 1 0,...
        1 0 0 1 1 0 0 1 1 1 0 0 1 0 0 0 0 1];

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
span = 10; % number of symbols for length of filter impulse response
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

phase_offset = pi/8;
transmitted_signal = modulate_carrier(transmited_baseband, Fc, t, phase_offset);

%% Transmit Through Wireless Channel

snr = 30;
received_signal = add_channel_impairments(transmitted_signal, B, Fc, Fs, snr);
% received_signal = transmitted_signal;

% Plot the Received Signal
figure;
plot(t, received_signal)
title("Received Signal")
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 0.0004]);

%% Demodulation

% Encorporate a frequency offset
frequency_offset = Fc * 0.1;
Fr = Fc + frequency_offset;

% Normalize the signal for stability of PID tuning constants
normalized_signal = normalize(received_signal);

% PID Tuning Constants
Kp = 0.1;
Ki = 0.002;
Kd = 0;

% Costas Loop to extract I and Q components
[I, Q, theta, err] = costas_loop(normalized_signal, Fc, Fs, t, Kp, Ki, Kd);

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

% TODO 1.1: Apply SPS frequency offset at the receiver
received_sps = ?

% TODO 1.2: Apply SPS phase offset at the receiver
received_samples = ?

% Naively sample the received signal at the symbol rate to get the received symbols
naive_received_symbols = downsample(received_samples, received_sps);

% Visualize the received symbols in a constellation diagram (scatterplot)
scatterplot(naive_received_symbols);
title("Naively Sampled Symbols");

% Normalize for stability for PID stability
received_samples = normalize(real(received_samples));

% TODO Section 2: Implement one (or multiple) of the timing error detectors to properly sample
% PID Tuning Constants
Kp = ? % observe how Kp causes faster settling w/o oscillation but increases self-noise effect
Ki = ? % observe how Ki smooths the oscillations before settling to a constant offset (but eventually goes unstable)
Kd = ? % observe how Kd increases high frequency noise in the lock

% Recreate the downsample function iteratively, then implement the TED

symbs = ?

% TODO: Visualize the properly receieved samples
% Note: A good chunk of the early symbols might be garbage as the TED locks on, so feel free
%       to only plot the later end of the symbols

% TODO: You'll probably want to change this once you have the timing error detector
received_symbols = naive_received_symbols;

% Apply Thresholding to detect the symbols (turning them back to 1's and -1's)
received_symbols(received_symbols > 0) = 1;
received_symbols(received_symbols <= 0) = -1;

%% Map Symbols back to Bits

% Map the detected symbols back to bits (turn them back to 1's and 0's
detected_bits = received_symbols;
detected_bits(detected_bits == -1) = 0;

%% Frame Syncronization!

key = [+1 +1 +1 -1 -1 -1 +1 -1 -1 +1 -1];

% Frame Sync Algorithm
num_message_chars = 7;

% TODO 3.1: Cross correlate and extract where the message starts

% TODO 3.2: Plot the cross correlation with the key


%% Extract Binary message

% TODO 3.3: Get binary message to decode
% message = ?

%% Decode and Display Message

% TODO 3.4: Uncomment this and you should see your final message!
% s = string(message);
% s = strjoin(s);
% s = strrep(s, " ", "");
% disp(s);
% s_len = length(char(s));
% inputString = char(s);
% binaryString = inputString(1:end-mod(length(inputString),8));
% binaryChunks = reshape(binaryString, 8, []).';
% asciiChars = char(bin2dec(binaryChunks)).';
% disp(asciiChars);
