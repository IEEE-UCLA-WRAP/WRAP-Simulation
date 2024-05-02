clear;
close all;

%%
r = 5; % Inverse of code rate

% Values that we're giving you
Fs = 4e6; % Sampling rate (samples/sec)
Rs = 50e3; % Symbol rate in symbols/sec (baud)
Tmax = 0.1 * r; % Max time for the simulation (sec)
fc = 1e6; % Carrier frequency (Hz)
B = 350e6; % if the 

% Complete these expressions using the variables above
N = Fs*Tmax; % Total number of sample points in the simulation
Ns = Rs*Tmax; % Number of symbols to send
sps = Fs/Rs; % Number of samples per symbol

% Use these variables for plotting
t = linspace(0,Tmax,N); % Time vector
f = linspace(-Fs/2, Fs/2, N); % Frequency vector. Recall from 113 that DFT gives us -f/2 to +f/2

%% Symbol Generation
% Generate random set of bits
original_bits = randi([0, 1], [1, Ns/r]);

% Repetion Code
bits = generate_repetition_code(original_bits,r);

% Map bits to BPSK symbols
symbs = 2 * (bits - 0.5);

scatterplot(symbs);

% upsample
samples = symbs;
deltas = upsample(samples,sps);

figure;
stem(deltas);
xlim([0 length(deltas)/100]);

%% Pulse Shaping (RRC)

% RRC setup
span = 5;           % number of symbols to truncate RRC
rolloff = 0.5;      % parameter of RRC
RRC = rcosdesign(rolloff, span, sps, 'sqrt');

rect_filt = ones([1, sps]);

symbol_wave = conv(deltas, rect_filt, 'same'); %%%%%%%%%%%%%%%%%%%%

figure;
stem(symbol_wave);
xlim([0 length(deltas)/100]);

%% Carrier Modulation
comp_carr = exp(-1i*2*pi*fc*t);
analytic_sig = symbol_wave .* comp_carr;

Tx = real(analytic_sig);
figure;
stem(Tx);
xlim([0 length(Tx)/1000]);

%% Noise and Band Limited Channel
% std_ = 0.01;
noise = 0;%std_ * randn(1,N);
Rx = Tx + noise;
% 
y = Rx;% apply_channel_bandlimit(Rx, fc, B, Fs); %%%%%%%%%%%%%%%%%%%%%%% Rx
%% SRRC

c = 2*cos(2*pi*fc*t);
s = -2*sin(2*pi*fc*t); %because we get I/2 and Q/2
% demodulate the signal by mixing with NCO signals
inph = y .* c;
quad = y .* s;
% low pass filter to get rid of high frequency components
inph = lowpass(inph, 1.5*10^6, Fs);%filter(lp, 1, inph); might need to LPF them when they're still together??
quad = lowpass(quad, 1.5*10^6, Fs);%filter(lp, 1, quad);

rec_bb = inph + 1i*quad;

y = upfirdn(rec_bb, RRC, 1, 1);
y = y(length(RRC)/2+0.5:end-length(RRC)/2+0.5);

rec_samples = y(1:sps:N);

figure();
scatterplot(rec_samples);
title("RX Constellation Diagram")

%% Bit conversion
%Conversion to bits
bits_received = double(rec_samples > 0);

% Decode bits
coded_bits_received = decode_repetition_code(bits_received, r);

bit_error_rate = compute_bit_error_rate(original_bits, coded_bits_received);
fprintf("Coded Bit error rate: %g %%\n",bit_error_rate);

bit_error_rate = compute_bit_error_rate(bits, bits_received);
fprintf("Bit error rate: %g %%\n",bit_error_rate);

%% Local functions
function bit_error_rate = compute_bit_error_rate(bits_intended, bits_received)
    if (length(bits_intended) ~= length(bits_received))
        bit_error_rate = -1;
    else
        bit_error_rate = (sum(bits_received ~= bits_intended)/length(bits_intended))*100;
    end
end

function received_signal = apply_channel_bandlimit(transmitted_signal, fc, B, Fs)
    s = tf('s');
    Hc = B*s/(s^2 + s*B + fc^2);
    Hd = c2d(Hc, 1/Fs, 'tustin');
    [a, b] = tfdata(Hd);
    
    received_signal = real(filter(a{:}, b{:}, transmitted_signal));
end

