%% Simulation Setup

Fc = 1e6; % Carrier frequency in Hz
Fs = 4e6; % Sampling rate
Rs = 50e3; % Symbol rate in symbols/sec (baud)
B = 40e3; % Bandwidth

sps = Fs / Rs; % Number of samples per symbol
num_symbols = 296; % Number of symbols/bits to send
num_samples = num_symbols * sps; % Number of discrete samples

% Use these guys for plotting
t = linspace(0, num_symbols / Rs, num_samples); % Time variable
f = linspace(-Fs/2,Fs/2, num_samples); % Frequency variable

% Clear all previous figures before starting
close all;

%% Mapping Bits to Symbols

% Below are the bits for the binary message that we want to send

bits = [0 1 0 1 0 1 1 1 0 1 1 0 1 1 1 1 0 1 1 1 0 1 1 1 0 0 1 0 0 0 0 0 0 1,... 
        0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 1 0 0 1 1 0 0 0 1 1 0 1 1 1 1 0 1 1 1,... 
        0 1 1 0 0 1 1 0 0 1 0 1 0 0 1 0 0 0 0 0 0 1 0 1 0 1 1 1 0 1 0 1 0 0,... 
        1 0 0 1 0 0 0 0 0 1 0 1 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1 0 0 1 1,... 
        0 1 1 0 1 1 1 1 0 0 1 0 0 0 0 0 0 1 0 0 1 1 0 1 0 1 1 1 0 1 0 1 0 1,... 
        1 0 0 0 1 1 0 1 1 0 1 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1 0 1 0 0 0 1 1 0,... 
        1 0 0 0 0 1 1 0 1 0 0 1 0 1 1 1 0 0 1 1 0 0 1 0 0 0 0 0 0 1 0 0 1 0,... 
        0 1 0 1 1 1 0 0 1 1 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1 1 0 1 1 0 1 1 1 1,... 
        0 1 1 0 1 1 1 1 0 1 1 0 1 1 0 0 0 0 1 0 0 0 0 1];

% PSA: For the outlined steps in the starter code, each step may take one or
% many lines of code to implement. I only included the "symbols = ?" in the next
% step (and similarly for future steps) so that variable names are consistent
% among students for easier debugging

% TODO 1.1.1: Convert to BPSK Symbols (1's and -1's)
symbols = bits;
symbols(symbols == 0) = -1;

% Visualize BPSK Symbols in a Constellation Diagram
scatterplot(symbols);
title("Transmitted Symbols");

%% Represent Symbols as Time Shifted Deltas

% TODO 1.2.1: Use upsample to form a delta train
deltas = upsample(symbols, sps);

% Visualize Deltas (xlimited to only first ~20 symbols)
figure;
stem(t, deltas);
title("Upsampled Deltas");
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 0.0004]);

%% Filter Deltas with Transmitter Filter

% TODO 1.3.1: Generate rectangular pulse shape filter
ps_filter = ones(1, sps);

% TODO 2.2: Replace the rect with an SRRC as the pulse shape filter 
% (ignore this if you are not at assignment section 2 yet)

% TODO 1.3.2: Convolve the deltas with the rectangular window
transmited_baseband = conv(deltas, ps_filter, 'same'); % convolve deltas with srrc

% Visualize Transmitted Baseband Signal
figure;
plot(t, transmited_baseband);
title("Transmited Baseband Signal");
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 0.0004]);

%% Modulate Baseband Signal to Passband

% TODO 1.4.1: Modulate the baseband signal to passband
complex_carrier = exp(-1i*2*pi*Fc*t);
analytic_signal = transmited_baseband .* complex_carrier;
transmitted_signal = real(analytic_signal);

% Plot the fft of the transmitted signal
figure;
plot(f, fftshift(abs(fft(transmitted_signal))));
title("FFT of Transmitted Signal");
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% TODO 1.4.2: Visualize Transmitted Passband Signal 
% (From now on, you do the plotting!)
figure;
plot(t, transmitted_signal)
title("Transmited Passband Signal")
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 0.0004]);

%% Transmit Through Wireless Channel

% TODO 2.1: Apply the channel bandlimit
% (ignore this if you are not at assignment section 2 yet)

% TODO 1.5.1: Add some AWGN to the signal
snr = 20;
received_signal = awgn(transmitted_signal, snr, 'measured');

% TODO 1.5.2: Plot the Corrupted Signal
figure;
plot(t, received_signal)
title("Received Corrupted Signal")
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 0.0004]);

%% Demodulation

% TODO 1.6.1: Start demodulation by mixing the received signal by mixing with the sinusoids
I_ = 2 * received_signal .* cos(2*pi*Fc*t);
Q_ = 2 * received_signal .* sin(2*pi*Fc*t);

% Lowpass filter the previous results to get I and Q
I = lowpass(I_, 1.5e6, Fs);
Q = lowpass(Q_, 1.5e6, Fs);

% TODO 1.6.2: Plot the In Phase Component
figure;
plot(t, I)
title("Received In Phase Component")
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 0.0004]);

%% Apply Receiver Impulse Response

% TODO 1.7.1: Combine I and Q components to construct the received baseband signal (I + jQ)
received_baseband = I + 1j * Q;

% TODO 1.7.2: Convolve the received baseband signal with the receiver filter
received_samples = conv(received_baseband, ps_filter, 'same');

%% Sample and Detect Symbols

% TODO 1.8.1: (Naively) sample the received signal at the symbol rate to get the received symbols
received_symbols = downsample(received_samples, sps);

% TODO 1.8.2: Visualize the received symbols in a constellation diagram (scatterplot)
scatterplot(received_symbols);
title("Received Symbols");

% TODO 1.8.3: Apply Thresholding to detect the symbols (turning them back to 1's and -1's)
received_symbols(received_symbols > 0) = 1;
received_symbols(received_symbols <= 0) = -1;

%% Map Symbols back to Bits

% TODO 1.9.1: Map the detected symbols back to bits (turn them back to 1's
% and 0's
detected_bits = received_symbols;
detected_bits(detected_bits == -1) = 0;

% Calculate Bit Error Rate
num_bit_errors = sum(detected_bits(1:length(bits)) ~= bits)
BER = 100 * num_bit_errors / length(bits)

%% Decode and Display Message

% If you have successfully decoded the message, you should see the message "Congrats! Simulation is complete"
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
